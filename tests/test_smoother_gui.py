import unittest
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from autoqy_core.epsilon import EpsilonResult
from autoqy_core.smoother import SpectralDataset
from autoqy_core.tools.smoother_gui import (
    _csv_filename,
    _export_csv_payload,
    _pack,
    _pack_epsilon,
    create_app,
)

try:
    from dash import html
    from dash._callback_context import context_value
    from dash._utils import AttributeDict
except ImportError:
    html = None


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
        opened = [component for component in details if component.open is True]
        self.assertEqual(len(opened), 1)
        self.assertIn("1 · Data", str(opened[0].to_plotly_json()))

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
            [], 1, [None], [1.0],
        )
        self.assertIsNone(preview_result[4])
        self.assertIsNotNone(preview_result[5])
        self.assertFalse(preview_result[6])

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
