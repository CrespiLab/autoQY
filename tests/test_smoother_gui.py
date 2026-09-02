import unittest
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from autoqy_core.epsilon import EpsilonResult
from autoqy_core.plot_style import ANALYSIS_TRACE_PALETTE
from autoqy_core.smoother import SpectralDataset
from autoqy_core.tools.smoother_gui import (
    _colors,
    _combine_loaded,
    _csv_filename,
    _export_csv_payload,
    _pack,
    _pack_epsilon,
    _spectrum_colors,
    _wavelength_slice,
    create_app,
)

try:
    from dash import dcc, html
    from dash._callback_context import context_value
    from dash._utils import AttributeDict
except ImportError:
    dcc = html = None


def _components(root):
    if isinstance(root, (list, tuple)):
        for child in root:
            yield from _components(child)
        return
    yield root
    children = getattr(root, "children", None)
    if children is not None:
        yield from _components(children)


class ProcessedAbsorbanceExportTests(unittest.TestCase):
    def setUp(self):
        self.dataset = SpectralDataset(
            np.array([400.0, 401.0, 402.0]),
            np.array([0.0]),
            np.array([[0.1], [0.2], [0.3]]),
            source_format="csv",
        )
        self.packed = _pack(self.dataset, ["irradiation"], ["irradiation.csv"])

    def test_processed_absorbance_has_a_csv_export_without_epsilon(self):
        filename, text = _export_csv_payload(None, self.packed)
        self.assertEqual(filename, "processed-absorbance.csv")
        self.assertTrue(text.startswith("Wavelength,0\n"))
        self.assertIn("2.00000000e-01", text)

    def test_negative_processed_absorbance_can_be_converted_to_zero(self):
        self.packed["absorbance"][1][0] = -0.2
        _, text = _export_csv_payload(None, self.packed, nonnegative=True)
        frame = pd.read_csv(StringIO(text))
        np.testing.assert_array_equal(frame["0"], [0.1, 0.0, 0.3])

    def test_negative_absorbance_and_epsilon_can_be_converted_to_zero(self):
        result = EpsilonResult(
            wavelengths=np.array([400.0, 401.0]),
            absorbance=np.array([[-0.1], [0.2]]),
            concentrations_m=np.array([1e-4]),
            path_lengths_cm=np.array([1.0]),
            individual=np.array([[-1000.0], [2000.0]]),
            mean=np.array([-1000.0, 2000.0]),
            standard_deviation=np.array([10.0, 20.0]),
            standard_error=np.array([10.0, 20.0]),
        )
        _, text = _export_csv_payload(
            _pack_epsilon(result, ["sample"]), None, nonnegative=True
        )
        frame = pd.read_csv(StringIO(text))
        for column in ("Absorbance__sample", "Epsilon_M-1_cm-1__sample",
                       "Epsilon_mean_M-1_cm-1"):
            self.assertGreaterEqual(frame[column].min(), 0.0)

    def test_custom_csv_name_accepts_an_optional_extension(self):
        self.assertEqual(_csv_filename("treated irradiation", "default.csv"),
                         "treated irradiation.csv")
        self.assertEqual(_csv_filename("treated.CSV", "default.csv"), "treated.CSV")
        self.assertEqual(_csv_filename("", "default.csv"), "default.csv")

    def test_custom_csv_name_rejects_paths_and_other_extensions(self):
        with self.assertRaisesRegex(ValueError, "without a folder path"):
            _csv_filename("subfolder/treated.csv", "default.csv")
        with self.assertRaisesRegex(ValueError, "must use the .csv extension"):
            _csv_filename("treated.tsv", "default.csv")


@unittest.skipUnless(html is not None, "Dash GUI dependencies are not installed")
class SpectralGuiTests(unittest.TestCase):
    def setUp(self):
        self.app = create_app()

    def test_only_data_panel_is_open_initially(self):
        details = []

        def visit(component):
            if isinstance(component, (list, tuple)):
                for child in component:
                    visit(child)
                return
            if isinstance(component, html.Details):
                details.append(component)
            children = getattr(component, "children", None)
            if children is not None:
                visit(children)

        visit(self.app.layout)
        opened = [
            component for component in details
            if getattr(component, "open", False) is True
        ]
        self.assertEqual(len(opened), 1)
        self.assertIn("1 · Data", str(opened[0].to_plotly_json()))

    def test_plot_controls_offer_legend_toggle_slice_and_exports(self):
        components = list(_components(self.app.layout))
        by_id = {
            component.id: component for component in components
            if getattr(component, "id", None)
        }
        self.assertNotIn("spectrum-legend-labels", by_id)
        self.assertEqual(by_id["show-spectrum-legend"].value, ["on"])
        self.assertEqual(by_id["minimal-spectrum-colors"].value, [])
        self.assertEqual(by_id["include-plot-header"].value, [])
        self.assertEqual(by_id["origin-epsilon-export"].value, [])
        self.assertEqual(by_id["origin-slice-export"].value, [])
        self.assertFalse(by_id["wavelength-slice-panel"].open)
        self.assertNotIn("responsive", by_id["epsilon-plot"].config)
        for component_id in (
            "save-epsilon-png", "save-epsilon-svg", "save-slice-png",
            "save-slice-svg", "save-slice-csv", "save-nmr-png", "save-nmr-svg",
        ):
            self.assertIn(component_id, by_id)
        self.assertIn("epsilon-image-message.children", self.app.callback_map)
        self.assertIn("nmr-image-message.children", self.app.callback_map)
        self.assertIn("window.showSaveFilePicker", self.app.index_string)
        self.assertIn("exportLayout.showlegend = false", self.app.index_string)
        self.assertIn("exportLayout.title = null", self.app.index_string)
        self.assertIn("width = 1200", self.app.index_string)
        self.assertIn("height = 900", self.app.index_string)

    def test_main_plot_uses_the_analysis_palette(self):
        self.assertEqual(_colors(), list(ANALYSIS_TRACE_PALETTE))

    def test_minimal_palette_marks_initial_intermediate_and_final_spectra(self):
        self.assertEqual(_spectrum_colors(1, True), ["#2d6f8e"])
        self.assertEqual(
            _spectrum_colors(4, True),
            [
                "#2d6f8e",
                "rgba(108,114,128,0.38)",
                "rgba(108,114,128,0.38)",
                "#d67b36",
            ],
        )
        self.assertEqual(_spectrum_colors(7), [
            ANALYSIS_TRACE_PALETTE[index % len(ANALYSIS_TRACE_PALETTE)]
            for index in range(7)
        ])

    def test_single_time_series_keeps_coordinates_and_interpolates_a_slice(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]),
            np.array([0.0, 2.5, 5.0]),
            np.array([[0.0, 1.0, 2.0], [10.0, 11.0, 12.0]]),
            source_format="spectragryph",
        )
        combined, _, _ = _combine_loaded([(dataset, "series.dat")])
        np.testing.assert_array_equal(combined.coordinates, dataset.coordinates)
        selected, coordinates, values = _wavelength_slice(combined, 405.0)
        self.assertEqual(selected, 405.0)
        np.testing.assert_array_equal(coordinates, dataset.coordinates)
        np.testing.assert_allclose(values, [5.0, 6.0, 7.0])

    def test_wavelength_slice_rejects_values_outside_the_processed_range(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0]),
            np.array([[0.0], [1.0]]), source_format="csv",
        )
        with self.assertRaisesRegex(ValueError, "between 400 and 410"):
            _wavelength_slice(dataset, 420.0)

    def test_slice_callback_builds_plot_export_data_and_custom_axes(self):
        callback = next(
            value["callback"].__wrapped__
            for key, value in self.app.callback_map.items()
            if "wavelength-slice-store.data" in key
        )
        dataset = SpectralDataset(
            np.array([400.0, 410.0]),
            np.array([0.0, 5.0]),
            np.array([[0.0, 2.0], [10.0, 12.0]]),
            source_format="spectragryph",
        )
        figure, data, message = callback(
            _pack(dataset, ["start", "end"], ["series.dat"]),
            405.0, "Elapsed time (s)", "Optical density",
        )
        self.assertEqual(data["values"], [5.0, 7.0])
        self.assertEqual(data["x_label"], "Elapsed time (s)")
        self.assertEqual(data["y_label"], "Optical density at 405 nm")
        self.assertEqual(figure.layout.xaxis.title.text, "Elapsed time (s)")
        self.assertEqual(figure.layout.yaxis.title.text, "Optical density")
        self.assertIn("405 nm", message)

    def test_preview_and_save_work_without_concentrations(self):
        callbacks = self.app.callback_map
        preview = next(
            value["callback"].__wrapped__
            for key, value in callbacks.items()
            if "processed-store.data" in key
        )
        save = next(
            value["callback"].__wrapped__
            for key, value in callbacks.items()
            if "epsilon-save-message" in key
        )
        dataset = SpectralDataset(
            np.array([400.0, 401.0, 402.0]),
            np.array([0.0]),
            np.array([[0.1], [-0.2], [0.3]]),
            source_format="csv",
        )
        packed = _pack(dataset, ["irradiation"], ["irradiation.csv"])
        preview_result = preview(
            packed, 400, 402, [], None, None, "off", 5, 3,
            [], 1, [None], [1.0], [], [], "Custom wavelength", "Custom OD",
            "Custom epsilon",
        )
        self.assertIsNone(preview_result[4])
        self.assertIsNotNone(preview_result[5])
        self.assertFalse(preview_result[6])
        self.assertFalse(preview_result[0].layout.showlegend)
        self.assertEqual(preview_result[0].data[0].name, "irradiation")
        self.assertEqual(preview_result[0].layout.xaxis.title.text, "Custom wavelength")
        self.assertEqual(preview_result[0].layout.yaxis.title.text, "Custom OD")
        self.assertEqual(preview_result[5]["labels"], ["irradiation"])

        with TemporaryDirectory() as temporary:
            token = context_value.set(AttributeDict(
                triggered_inputs=[{"prop_id": "export-epsilon.n_clicks"}]
            ))
            try:
                save_result = save(
                    1, None, None, preview_result[5], temporary,
                    "my-treated-irradiation", ["on"],
                )
            finally:
                context_value.reset(token)
            destination = Path(temporary) / "my-treated-irradiation.csv"
            self.assertTrue(destination.is_file())
            self.assertIn(str(destination), save_result[0])
            frame = pd.read_csv(destination)
            self.assertEqual(frame["0"].min(), 0.0)


if __name__ == "__main__":
    unittest.main()
