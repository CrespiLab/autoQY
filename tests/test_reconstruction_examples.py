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
from autoqy_core.smoother import load_spectral_bytes
from autoqy_core.tools.smoother_gui import (
    _csv_exports, _load_local_paths, _load_reconstruction_epsilon, _pack_epsilon,
    _prepare_processing, _reconstruct_epsilon_pair, create_app,
    _interpolate_to_axis, _process_nmr_pair,
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
        self.raw_pair = {"filename": "reactant_pss.dat",
                         "reactant_wavelengths": raw.wavelengths.tolist(),
                         "reactant_absorbance": raw.absorbance[:, 0].tolist(),
                         "pss_wavelengths": raw.wavelengths.tolist(),
                         "pss_absorbance": raw.absorbance[:, -1].tolist()}
        self.calculate = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                              if "nmr-result-store.data" in key)

    def test_original_395nm_calibration_matches_bundled_product(self):
        request = {"mode": "calibrated", "percent": 77, "error_percent": 1, "error_basis": "sd"}
        result = self.calculate(self.raw_pair, self.epsilon_packed, 77, 1,
                                ["on"], 600, 650, "savgol", 5, 3, "calibrated", "sd", request)
        self.assertIsNotNone(result[3], result[2])
        target = self.epsilon.wavelengths
        reference = _interpolate_to_axis(self.raw_pair["reactant_wavelengths"], self.raw_pair["reactant_absorbance"], target)
        mixture = _interpolate_to_axis(self.raw_pair["pss_wavelengths"], self.raw_pair["pss_absorbance"], target)
        reference, mixture, _ = _process_nmr_pair(target, reference, mixture, ["on"], 600, 650, "savgol", 5, 3)
        legacy = reconstruct_product_from_nmr(target, reference, mixture, self.epsilon, 23, 1)
        for field in ("product", "product_lower", "product_upper"):
            np.testing.assert_array_equal(result[3][field], getattr(legacy, field))
        expected = pd.read_csv(EXAMPLE / "generic_inputs" / "product_absorptivity.csv")
        np.testing.assert_allclose(result[3]["wavelengths"], expected["Wavelength_nm"], rtol=1e-7)
        for field, column in (("product", "Product_epsilon_raw_M-1_cm-1"),
                              ("product_lower", "Product_lower_nonnegative_M-1_cm-1"),
                              ("product_upper", "Product_upper_nonnegative_M-1_cm-1")):
            # The bundled precomputed artifact differs slightly from the current
            # unmodified calibration pipeline, especially near a zero baseline.
            # Assert the absolute difference, not a relative error near zero.
            self.assertLess(np.max(np.abs(np.asarray(result[3][field]) - expected[column])), .8)

    def epsilon_pair(self):
        # A mathematical self-mixture sanity check using the real example CSV.
        # This is not a second independent PSS measurement or a new product estimate.
        text = (EXAMPLE / "generic_inputs" / "reactant_absorptivity.csv").read_text(encoding="utf-8")
        content = "data:text/csv;base64," + base64.b64encode(text.encode()).decode()
        source = _load_reconstruction_epsilon(content, "reactant_absorptivity.csv")
        return {"kind": "epsilon", "reference": source, "mixture": deepcopy(source)}

    def test_reference_triplicate_uv_pss_nmr_without_pss_errors(self):
        reference = self.epsilon_pair()["reference"]
        pair = self.raw_pair | {"kind": "calibrated", "reference": reference}
        request = {"mode": "calibrated", "percent": 77, "error_percent": 1, "error_basis": "sd"}
        sd = self.calculate(pair, None, 77, 1, ["on"], 600, 650, "savgol", 5, 3,
                            "calibrated", "sd", request)
        self.assertIsNotNone(sd[3], sd[2])
        self.assertIn("PSS error excluded", sd[3]["uncertainty_method"])
        sem = self.calculate(pair, None, 77, 1, ["on"], 600, 650, "savgol", 5, 3,
                             "calibrated", "sem", request | {"error_basis": "sem"})
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
        snapshot = {"version": 1, "dataset": None, "nmr": pair, "isomer_request": request,
                    "settings": {"reconstruction-mode": "calibrated", "isomer-percent": 77,
                                 "nmr-error-percent": 1, "reconstruction-error-basis": "sd"}}
        self.assertEqual(decode_session(encode_session(snapshot)), snapshot)

    def test_example_epsilon_upload_standalone_reconstruction_export_and_recovery(self):
        pair = self.epsilon_pair()
        request = {"mode": "epsilon", "percent": 100, "error_percent": 0, "error_basis": "sd"}
        result = self.calculate(pair, None, 100, 0, [], None, None, "off", 5, 3,
                                "epsilon", "sd", request)
        self.assertIsNotNone(result[3], result[2])
        data = result[3]
        np.testing.assert_allclose(data["product"], pair["mixture"]["epsilon"]["mean"])
        np.testing.assert_allclose(data["product_error"], pair["mixture"]["epsilon"]["sd"])
        ready = self.app.callback_map["export-epsilon.disabled"]["callback"].__wrapped__
        self.assertFalse(ready("isomer", None, None, data, request, None))
        name, text = _csv_exports("isomer", None, None, isomer=request, nmr_data=data)[0]
        destination = Path(self.temp.name) / name
        destination.write_text(text, encoding="utf-8")
        envelope = load_epsilon_envelope(destination)
        np.testing.assert_allclose(envelope.upper, data["product_upper"], rtol=1e-7)
        for basis in ("sd", "sem"):
            product = _reconstruct_epsilon_pair(pair, 77, 1, basis)
            stable, mixture = np.asarray(product["reference"]), np.asarray(product["mixture"])
            a = np.asarray(product["reference_error"])
            b = np.asarray(product["mixture_error"])
            expected = np.sqrt((b / .77)**2 + (.23 * a / .77)**2 +
                               ((stable - mixture) * .01 / .77**2)**2)
            np.testing.assert_allclose(product["product_error"], expected)
        snapshot = {"version": 1, "dataset": None, "nmr": pair,
                    "isomer_request": request, "settings": {"reconstruction-mode": "epsilon"}}
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
