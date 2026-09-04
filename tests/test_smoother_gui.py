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
    _append_packed,
    _absorbance_figure,
    _colors,
    _combine_loaded,
    _csv_filename,
    _export_csv_payload,
    _loaded_spectrum_rows,
    _pack,
    _pack_epsilon,
    _remove_packed,
    _reorder_packed,
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
        self.assertNotIn("show-spectrum-legend", by_id)
        self.assertIn("loaded-spectrum-manager", by_id)
        self.assertIn("show-all-legends", by_id)
        self.assertIn("hide-all-legends", by_id)
        self.assertEqual(by_id["minimal-spectrum-colors"].value, [])
        self.assertEqual(by_id["include-plot-header"].value, [])
        self.assertEqual(by_id["origin-epsilon-export"].value, [])
        self.assertEqual(by_id["origin-slice-export"].value, [])
        self.assertEqual(by_id["show-epsilon-export-grid"].value, [])
        self.assertEqual(by_id["show-slice-export-grid"].value, [])
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
        self.assertIn("width = 1300", self.app.index_string)
        self.assertIn("height = 1000", self.app.index_string)
        self.assertIn("axis.ticks = 'outside'", self.app.index_string)
        self.assertIn("axis.mirror = true", self.app.index_string)
        self.assertIn("axis.linewidth = 3", self.app.index_string)
        self.assertIn("Math.max(3.5, lineWidth * 1.8)", self.app.index_string)
        self.assertIn("size: 32", self.app.index_string)
        self.assertIn("showgrid: includeGrid", self.app.index_string)

    def test_main_plot_uses_the_analysis_palette(self):
        self.assertEqual(_colors(), list(ANALYSIS_TRACE_PALETTE))

    def test_legend_can_show_only_selected_spectrum_traces(self):
        import plotly.graph_objects as go

        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 1.0, 2.0]),
            np.array([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]]),
            source_format="csv",
        )
        figure = _absorbance_figure(
            go, dataset, dataset.absorbance, dataset.absorbance,
            ["one", "two", "three"], "off",
            legend_visibility=[True, False, True],
        )
        self.assertEqual(
            [trace.showlegend for trace in figure.data], [True, False, True]
        )
        self.assertTrue(figure.layout.showlegend)

    def test_loaded_spectrum_rows_include_editable_legend_names(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 1.0]),
            np.array([[1.0, 2.0], [1.5, 2.5]]), source_format="csv",
        )
        rows = _loaded_spectrum_rows(
            html, dcc, _pack(dataset, ["one", "two"], ["one.csv", "two.csv"])
        )
        components = list(_components(rows))
        legend_inputs = [
            component for component in components
            if isinstance(getattr(component, "id", None), dict)
            and component.id.get("type") == "legend-name"
        ]
        self.assertEqual([component.value for component in legend_inputs], ["one", "two"])

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

    def test_multiple_files_are_naturally_sorted_by_filename(self):
        def spectrum(value):
            return SpectralDataset(
                np.array([400.0, 410.0]), np.array([0.0]),
                np.array([[value], [value + 0.5]]), source_format="csv",
            )

        combined, labels, _ = _combine_loaded([
            (spectrum(10.0), "sample10.csv"),
            (spectrum(2.0), "sample2.csv"),
            (spectrum(1.0), "Sample1.csv"),
            (spectrum(3.0), "sampleA.csv"),
            (spectrum(4.0), "sampleB.csv"),
        ])
        self.assertEqual(
            labels, ["Sample1", "sample2", "sample10", "sampleA", "sampleB"]
        )
        np.testing.assert_array_equal(
            combined.absorbance[0], [1.0, 2.0, 10.0, 3.0, 4.0]
        )

    def test_repeated_single_spectrum_drops_append_in_elapsed_time(self):
        first = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0]),
            np.array([[1.0], [2.0]]), source_format="csv",
        )
        second = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0]),
            np.array([[3.0], [4.0]]), source_format="csv",
        )
        combined = _append_packed(
            _pack(first, ["spectrum"], ["spectrum.csv"]),
            _pack(second, ["spectrum"], ["spectrum.csv"]),
            elapsed_seconds=2.5,
        )
        self.assertEqual(combined["labels"], ["spectrum", "spectrum (2)"])
        self.assertEqual(combined["coordinates"], [0.0, 2.5])
        np.testing.assert_array_equal(
            np.asarray(combined["absorbance"]),
            np.array([[1.0, 3.0], [2.0, 4.0]]),
        )

    def test_selected_spectra_can_be_removed_without_retiming_the_rest(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 2.5, 8.0]),
            np.array([[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]]),
            source_format="csv",
        )
        packed = _pack(
            dataset, ["first", "middle", "last"], ["series.csv"],
            concentrations=[1.0, 2.0, 3.0], path_lengths=[1.0, 1.0, 2.0],
        )
        reduced = _remove_packed(
            packed, [1], concentrations=[10.0, 20.0, 30.0],
            path_lengths=[0.5, 1.0, 1.5],
            legend_values=[["on"], [], ["on"]],
            legend_names=["Start", "Middle", "Finish"],
        )
        self.assertEqual(reduced["labels"], ["first", "last"])
        self.assertEqual(reduced["coordinates"], [0.0, 8.0])
        self.assertEqual(reduced["concentrations"], [10.0, 30.0])
        self.assertEqual(reduced["path_lengths"], [0.5, 1.5])
        self.assertEqual(reduced["legend_visibility"], [True, True])
        self.assertEqual(reduced["legend_names"], ["Start", "Finish"])
        np.testing.assert_array_equal(
            np.asarray(reduced["absorbance"]),
            np.array([[1.0, 5.0], [2.0, 6.0]]),
        )
        self.assertIsNone(_remove_packed(packed, [0, 1, 2]))

    def test_loaded_spectra_can_be_reordered_with_their_settings(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 4.0, 9.0]),
            np.array([[1.0, 2.0, 3.0], [11.0, 12.0, 13.0]]),
            source_format="csv",
        )
        packed = _pack(
            dataset, ["one", "two", "three"],
            ["one.csv", "two.csv", "three.csv"],
        )
        reordered = _reorder_packed(
            packed, [1, 0, 2], concentrations=[1.0, 2.0, 3.0],
            path_lengths=[0.5, 1.0, 1.5],
            legend_values=[["on"], [], ["on"]],
            legend_names=["Initial", "Intermediate", "Final"],
        )
        self.assertEqual(reordered["labels"], ["two", "one", "three"])
        self.assertEqual(reordered["filenames"], ["two.csv", "one.csv", "three.csv"])
        self.assertEqual(reordered["coordinates"], [0.0, 4.0, 9.0])
        self.assertEqual(reordered["concentrations"], [2.0, 1.0, 3.0])
        self.assertEqual(reordered["legend_visibility"], [False, True, True])
        self.assertEqual(
            reordered["legend_names"], ["Intermediate", "Initial", "Final"]
        )
        np.testing.assert_array_equal(
            np.asarray(reordered["absorbance"])[0], [2.0, 1.0, 3.0]
        )

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
            [], 1, [None], [1.0], [[]], ["Custom legend"], [],
            "Custom wavelength", "Custom OD", "Custom epsilon",
        )
        self.assertIsNone(preview_result[4])
        self.assertIsNotNone(preview_result[5])
        self.assertFalse(preview_result[6])
        self.assertFalse(preview_result[0].layout.showlegend)
        self.assertEqual(preview_result[0].data[0].name, "Custom legend")
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
