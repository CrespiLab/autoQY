"""Merged reconstruction regressions using the bundled 395 nm source data."""
import base64
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd

from autoqy_core.epsilon import (
    calculate_epsilon_statistics, export_epsilon_csv, load_epsilon_table,
    reconstruct_spectrum_with_uncertainty, reconstruct_product_from_nmr,
)
from autoqy_core.smoother import load_spectral_bytes, SpectralDataset
from autoqy_core.tools.smoother_gui import (
    _csv_exports, _load_local_paths, _load_reconstruction_epsilon, _pack_epsilon,
    _prepare_processing, _reconstruct_epsilon_pair, create_app,
    _interpolate_to_axis, _subset_epsilon, _pack, _unpack, _reconstruct_from_main,
)
from autoqy_core.tools.spectral_recovery import encode_session, decode_session
from autoqy_core.epsilon_uncertainty import load_epsilon_envelope

EXAMPLE = Path(__file__).parents[1] / "ExampleData" / "Example-4_395nm-EpsilonError"


class ReconstructionExamplesTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.app = create_app(Path(self.temp.name) / "autosave.aqys")
        source = EXAMPLE / "crespi_group_inputs" / "source_data"
        packed, _ = _load_local_paths(sorted(source.glob("*.Abs8")))
        dataset, _, values, _ = _prepare_processing(
            packed, 250, 700, ["on"], 600, 650, "savgol", 5, 3)
        self.epsilon = calculate_epsilon_statistics(
            dataset.wavelengths, values, [6.82e-5, 6.87e-5, 6.80e-5], [1, 1, 1])
        self.epsilon_packed = _pack_epsilon(self.epsilon, packed["labels"])
        raw = load_spectral_bytes((source / "reactant_pss.dat").read_bytes(), "auto")
        packed_uv = _pack(raw, [f"UV {i}" for i in range(raw.absorbance.shape[1])], ["reactant_pss.dat"])
        treated, _, uv_values, _ = _prepare_processing(
            packed_uv, 250, 700, ["on"], 600, 650, "savgol", 5, 3)
        self.main = _pack(SpectralDataset(treated.wavelengths, treated.coordinates,
                                         uv_values, treated.source_format, 0),
                          packed_uv["labels"], packed_uv["filenames"])
        self.request = {"reference": 0, "mixture": len(packed_uv["labels"]) - 1,
                        "mode": "calibrated", "percent": 77, "error_percent": 1, "error_basis": "sd"}
        self.raw_pair = {"filename": "reactant_pss.dat",
                         "reactant_wavelengths": raw.wavelengths.tolist(),
                         "reactant_absorbance": raw.absorbance[:, 0].tolist(),
                         "pss_wavelengths": raw.wavelengths.tolist(),
                         "pss_absorbance": raw.absorbance[:, -1].tolist()}
        self.calculate = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                              if "nmr-result-store.data" in key)

    def test_original_395nm_calibration_matches_bundled_product(self):
        result = self.calculate(self.main, {"reference": {"epsilon": self.epsilon_packed}}, self.request, "sd")
        self.assertIsNotNone(result[3], result[2])
        treated = _unpack(self.main)
        mask = ((self.epsilon.wavelengths >= treated.wavelengths.min()) &
                (self.epsilon.wavelengths <= treated.wavelengths.max()))
        epsilon = _subset_epsilon(self.epsilon, mask)
        target = epsilon.wavelengths
        reference = _interpolate_to_axis(treated.wavelengths, treated.absorbance[:, 0], target)
        mixture = _interpolate_to_axis(treated.wavelengths, treated.absorbance[:, -1], target)
        legacy = reconstruct_product_from_nmr(target, reference, mixture, epsilon, 23, 1)
        for field in ("product", "product_lower", "product_upper"):
            np.testing.assert_array_equal(result[3][field], getattr(legacy, field))
        expected = pd.read_csv(EXAMPLE / "generic_inputs" / "product_absorptivity.csv")
        expected = expected.iloc[np.flatnonzero(mask)]
        np.testing.assert_allclose(result[3]["wavelengths"], expected["Wavelength_nm"], rtol=1e-7)
        for field, column in (("product", "Product_epsilon_raw_M-1_cm-1"),
                              ("product_lower", "Product_lower_nonnegative_M-1_cm-1"),
                              ("product_upper", "Product_upper_nonnegative_M-1_cm-1")):
            # Main preprocessing now happens on the native UV grid, BEFORE
            # interpolation, unlike the bundled legacy artifact. Require agreement
            # within 0.1% of peak ε, without unstable relative errors near zero.
            self.assertLess(np.max(np.abs(np.asarray(result[3][field]) - expected[column])),
                            .001 * np.max(np.abs(expected[column])))

    def test_main_selection_is_used_and_reference_does_not_replace_absorbance(self):
        main = deepcopy(self.main)
        values = np.asarray(main["absorbance"])
        main["absorbance"] = np.column_stack([values[:, -1], values[:, 0], values[:, -1]]).tolist()
        main["coordinates"] = [0, 1, 2]
        main["labels"] = ["other", "initial", "PSS"]
        main["legend_names"] = ["Other curve", "Edited stable", "Edited PSS"]
        calibration = {"reference": {"epsilon": self.epsilon_packed}}
        selected = _reconstruct_from_main(main, calibration, self.request | {"reference": 1, "mixture": 2})
        original = _reconstruct_from_main(self.main, calibration, self.request)
        np.testing.assert_array_equal(selected["product"], original["product"])
        self.assertEqual(selected["stable_label"], "Edited stable")
        self.assertEqual(selected["pss_label"], "Edited PSS")
        with self.assertRaisesRegex(ValueError, "different spectra"):
            _reconstruct_from_main(main, calibration, self.request | {"reference": 1, "mixture": 1})

    def test_reconstruction_plot_only_has_epsilon_and_respects_minimal_colors(self):
        calibration = {"reference": {"epsilon": self.epsilon_packed}}
        result = self.calculate(self.main, calibration, self.request, "sem", ["on"], ["on"])
        self.assertIsNotNone(result[3], result[2])
        figure = result[0]
        traces = {trace.name: trace for trace in figure.data if trace.name}
        self.assertNotIn("yaxis2", figure.layout.to_plotly_json())
        self.assertFalse(any("normalized" in (trace.name or "").lower() or
                             "processed" in (trace.name or "").lower() for trace in figure.data))
        self.assertEqual(traces["Stable ε"].line.color, "#2d6f8e")
        self.assertEqual(traces["Product ε"].line.color, "#d67b36")
        self.assertEqual(traces["PSS reconstructed ε"].line.color, "#000000")
        for name in ("Stable ε", "PSS reconstructed ε", "Product ε"):
            self.assertEqual(traces[name].line.dash, "solid")
        self.assertTrue(all(shape.line.dash == "solid" for shape in figure.layout.shapes))
        self.assertTrue(result[2].startswith("Metastable ε ready."))
        self.assertNotIn("calibration", result[2])
        self.assertEqual(result[6], "")
        self.assertIn("Stable ε SEM envelope", traces)
        self.assertFalse(any(trace.line.color == "#8a5a9b" for trace in figure.data))
        hidden = self.calculate(self.main, calibration, self.request, "sd", [], ["on"])[0]
        self.assertTrue(all(trace.visible == "legendonly" for trace in hidden.data
                            if trace.name == "Product ε" or trace.legendgroup == "product-error"))

    def test_reference_errors_and_overlap_are_validated_for_main_workflow(self):
        calibration = {"reference": {"epsilon": deepcopy(self.epsilon_packed)}}
        calibration["reference"]["epsilon"]["sd"][0] = None
        with self.assertRaisesRegex(ValueError, "finite, non-negative SD"):
            _reconstruct_from_main(self.main, calibration, self.request)
        calibration["reference"]["epsilon"] = deepcopy(self.epsilon_packed)
        calibration["reference"]["epsilon"]["wavelengths"] = (
            np.asarray(self.epsilon_packed["wavelengths"]) + 1000).tolist()
        with self.assertRaisesRegex(ValueError, "overlapping wavelengths"):
            _reconstruct_from_main(self.main, calibration, self.request)

    def test_simplified_layout_and_scrollable_sheet(self):
        def components(node):
            if isinstance(node, (list, tuple)):
                for child in node:
                    yield from components(child)
            elif hasattr(node, "to_plotly_json"):
                yield node
                yield from components(getattr(node, "children", None))
        by_id = {item.id: item for item in components(self.app.layout)
                 if isinstance(getattr(item, "id", None), str)}
        self.assertNotIn("nmr-upload", by_id)
        self.assertNotIn("nmr-baseline-enabled", by_id)
        self.assertNotIn("reconstruction-mode", by_id)
        self.assertIn("isomer-reference", by_id)
        self.assertIn("epsilon-reference-upload", by_id)
        table = by_id["displayed-data-sheet"]
        self.assertEqual(table.page_action, "none")
        self.assertTrue(table.virtualization)
        self.assertEqual(table.fixed_rows, {"headers": True})

    def epsilon_pair(self):
        # A mathematical self-mixture sanity check using the real example CSV.
        # This is not a second independent PSS measurement or a new product estimate.
        text = (EXAMPLE / "generic_inputs" / "reactant_absorptivity.csv").read_text(encoding="utf-8")
        content = "data:text/csv;base64," + base64.b64encode(text.encode()).decode()
        source = _load_reconstruction_epsilon(content, "reactant_absorptivity.csv")
        return {"kind": "epsilon", "reference": source, "mixture": deepcopy(source)}

    def test_reference_triplicate_uv_pss_nmr_without_pss_errors(self):
        reference = self.epsilon_pair()["reference"]
        pair = {"reference": reference}
        request = self.request
        sd = self.calculate(self.main, pair, request, "sd")
        self.assertIsNotNone(sd[3], sd[2])
        self.assertIn("PSS error excluded", sd[3]["uncertainty_method"])
        sem = self.calculate(self.main, pair, request | {"error_basis": "sem"}, "sem")
        self.assertIsNotNone(sem[3], sem[2])
        np.testing.assert_array_equal(sd[3]["product"], sem[3]["product"])
        self.assertTrue(np.all(np.asarray(sem[3]["product_upper"]) <= np.asarray(sd[3]["product_upper"]) + 1e-10))
        self.assertTrue(np.all(np.asarray(sem[3]["product_lower"]) >= np.asarray(sd[3]["product_lower"]) - 1e-10))
        for kind in ("isomer", "nmr"):
            exports = _csv_exports(kind, None, None, isomer=request, nmr_data=sd[3])
            self.assertEqual(len(exports), 1 if kind == "isomer" else 2)
            product_name, text = exports[-1]
            destination = Path(self.temp.name) / product_name
            destination.write_text(text, encoding="utf-8")
            envelope = load_epsilon_envelope(destination)
            np.testing.assert_allclose(envelope.upper, sd[3]["product_upper"], rtol=1e-7)
        snapshot = {"version": 1, "dataset": self.main, "nmr": pair, "isomer_request": request,
                    "settings": {"isomer-percent": 77,
                                 "nmr-error-percent": 1, "reconstruction-error-basis": "sd"}}
        self.assertEqual(decode_session(encode_session(snapshot)), snapshot)

    def test_example_epsilon_upload_requires_main_spectra_and_zero_nmr_error_is_respected(self):
        pair = {"reference": self.epsilon_pair()["reference"]}
        request = self.request | {"error_percent": 0}
        absent = self.calculate(None, pair, request)
        self.assertIsNone(absent[3])
        self.assertIn("main window", absent[2])
        original = deepcopy(self.main)
        data = _reconstruct_from_main(self.main, pair, request)
        self.assertEqual(self.main, original)
        self.assertEqual(data["composition_error"], 0)
        with_nmr = _reconstruct_from_main(self.main, pair, self.request)
        self.assertTrue(np.all(np.asarray(data["product_upper"]) <= np.asarray(with_nmr["product_upper"]) + 1e-10))
        ready = self.app.callback_map["export-epsilon.disabled"]["callback"].__wrapped__
        self.assertFalse(ready("isomer", None, None, data, request, None))
        name, text = _csv_exports("isomer", None, None, isomer=request, nmr_data=data)[0]
        destination = Path(self.temp.name) / name
        destination.write_text(text, encoding="utf-8")
        envelope = load_epsilon_envelope(destination)
        np.testing.assert_allclose(envelope.upper, data["product_upper"], rtol=1e-7)
        snapshot = {"version": 1, "dataset": self.main, "nmr": pair,
                    "isomer_request": request, "settings": {"nmr-error-percent": 0}}
        self.assertEqual(decode_session(encode_session(snapshot)), snapshot)

    def test_missing_errors_and_nonoverlap_are_not_silently_accepted(self):
        pair = self.epsilon_pair()
        pair["mixture"]["epsilon"]["sd"][0] = None
        with self.assertRaisesRegex(ValueError, "unknown errors"):
            _reconstruct_epsilon_pair(pair, 77, 1)
        pair = self.epsilon_pair()
        pair["mixture"]["epsilon"]["wavelengths"] = (np.asarray(
            pair["mixture"]["epsilon"]["wavelengths"]) + 1000).tolist()
        with self.assertRaisesRegex(ValueError, "in common"):
            _reconstruct_epsilon_pair(pair, 77, 1)

    def test_independent_error_formula_and_monte_carlo_check(self):
        value, error = reconstruct_spectrum_with_uncertainty([10], [20], 77, [1], [2], 1)
        expected = np.sqrt((2/.77)**2 + (.23/.77)**2 + (10*.01/.77**2)**2)
        self.assertAlmostEqual(error[0], expected)
        rng = np.random.default_rng(3)
        reference, mixture, fraction = rng.normal(10, 1, 100000), rng.normal(20, 2, 100000), rng.normal(.77, .01, 100000)
        samples = (mixture - (1-fraction)*reference)/fraction
        self.assertLess(abs(samples.std() / expected - 1), .02)
        self.assertLess(abs(samples.mean() / value[0] - 1), .01)
