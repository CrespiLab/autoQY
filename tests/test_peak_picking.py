import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import plotly.graph_objects as go

from autoqy_core.smoother import SpectralDataset
from autoqy_core.tools.smoother_gui import (
    _figure_sheet_page, _overlay_peaks, _pack, _wavelength_peaks, create_app,
)
from autoqy_core.tools.spectral_recovery import SpectralRecovery


class PeakPickingTests(unittest.TestCase):
    def test_prominence_filters_weaker_maximum(self):
        x = np.arange(7.)
        y = [0, 10, 0, 1, 0, 0, 0]
        self.assertEqual(_wavelength_peaks(x, y, 5, 0).tolist(), [1, 3])
        self.assertEqual(_wavelength_peaks(x, y, 20, 0).tolist(), [1])

    def test_distance_uses_nm_on_irregular_descending_axis(self):
        x = np.array([390, 400, 400.2, 400.4, 405, 410, 420.])
        y = np.array([0, 10, 0, 8, 0, 7, 0.])
        self.assertEqual(_wavelength_peaks(x, y, 0, 1).tolist(), [1, 5])
        self.assertEqual(_wavelength_peaks(x[::-1], y[::-1], 0, 1).tolist(), [5, 1])

    def test_plateau_flat_short_endpoints_and_gaps(self):
        self.assertEqual(_wavelength_peaks(range(6), [0, 2, 2, 2, 0, 0], 0, 0).tolist(), [2])
        for y in ([1, 1, 1], [1, 2], [1, 2, 3], [3, 2, 1],
                  [0, 2, np.nan, 1, 0]):
            self.assertEqual(_wavelength_peaks(range(len(y)), y, 0, 0).tolist(), [])

    def test_invalid_controls(self):
        for prominence, separation in ((None, 1), (-1, 1), (101, 1),
                                       (np.nan, 1), (5, -1), (5, np.inf)):
            with self.assertRaises((ValueError, TypeError)):
                _wavelength_peaks([1, 2, 3], [0, 1, 0], prominence, separation)

    def test_overlay_preserves_table_and_respects_axes_and_visibility(self):
        figure = go.Figure()
        figure.add_trace(go.Scatter(x=[400, 401, 402], y=[0, 1, 0], mode="lines",
                                   name="Sample", yaxis="y2", line={"color": "blue"}))
        figure.add_trace(go.Scatter(x=[400, 401, 402], y=[0, 2, 0], mode="lines",
                                   name="Hidden", visible="legendonly"))
        figure.add_trace(go.Scatter(x=[400, 401, 402], y=[0, 3, 0], mode="lines",
                                   name="Error bound", line={"width": 0}))
        before = _figure_sheet_page(figure, size=None)
        _overlay_peaks(go, figure, ["on"], 5, 5)
        self.assertEqual(before, _figure_sheet_page(figure, size=None))
        self.assertEqual(len(figure.data), 4)
        marker = figure.data[-1]
        self.assertEqual(marker.text, ("401 nm",))
        self.assertEqual(marker.yaxis, "y2")
        self.assertEqual(marker.legendgroup, figure.data[0].legendgroup)
        self.assertIn("1 wavelength maxima", figure.layout.meta["peak_message"])

    def test_preview_toggle_recovery_and_exports(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / "session.aqys"
            SpectralRecovery(path).write({"version": 1, "settings": {
                "peak-picking-enabled": ["on"], "peak-prominence": 12,
                "peak-separation": 2}})
            app = create_app(path)
            def components(root):
                if isinstance(root, (list, tuple)):
                    for child in root:
                        yield from components(child)
                else:
                    yield root
                    yield from components(getattr(root, "children", []) or [])
            controls = {getattr(c, "id", None): c for c in components(app.layout)
                        if isinstance(getattr(c, "id", None), str)}
            self.assertEqual(controls["peak-prominence"].value, 12)
            self.assertEqual(controls["peak-picking-enabled"].value, ["on"])
            popup = controls["peak-picking-panel"].children[0].children[1]
            self.assertEqual(popup.className, "info-popup")
            range_panel = next(c for c in components(app.layout)
                               if getattr(c, "className", None) == "panel tool-details"
                               and "wavelength-low" in {getattr(child, "id", None)
                                                        for child in components(c)
                                                        if isinstance(getattr(child, "id", None), str)})
            self.assertIn("peak-picking-enabled", [getattr(c, "id", None)
                                                   for c in components(range_panel)])
            displayed_panel = next(c for c in components(app.layout)
                                   if getattr(c, "className", None) == "panel result-summary")
            self.assertNotIn("peak-picking-enabled", [getattr(c, "id", None)
                                                       for c in components(displayed_panel)])
            self.assertFalse(any("peak" in str(o) for o in controls["csv-export-kind"].options))
            callback = next(v["callback"].__wrapped__ for k, v in app.callback_map.items()
                            if "epsilon-plot.figure" in k)
            data = _pack(SpectralDataset(np.array([400., 401., 402.]), np.array([0.]),
                                        np.array([[0.], [1.], [0.]]), "csv"), ["Sample"], ["sample.csv"])
            args = (data, 400, 402, [], None, None, "off", 5, 3, [], None,
                    [], [], [["on"]], ["Sample"], [], [], None, None, None,
                    [{"type": "legend-spectrum", "index": 0}],
                    [{"type": "legend-name", "index": 0}], None, [])
            off = callback(*args)
            on = callback(*args, ["on"], 5, 5)
            self.assertEqual(on[-1], "")
            self.assertEqual(len(on[0].data), len(off[0].data) + 1)
            self.assertEqual(off[5], on[5])
            self.assertEqual(_figure_sheet_page(off[0], size=None),
                             _figure_sheet_page(on[0], size=None))
            self.assertEqual(len(_figure_sheet_page(on[0], size=None)[1]), 3)
            invalid = callback(*args, ["on"], None, 5)
            self.assertEqual(invalid[-1], "")
            self.assertIn("Peak picking:", invalid[0].layout.meta["peak_message"])
            epsilon_args = list(args)
            epsilon_args[11:13] = [[0.001], [1.]]
            epsilon = callback(*epsilon_args, ["on"], 5, 5)
            self.assertEqual(epsilon[-1], "")
            self.assertTrue(any(t.meta and t.meta.get("peak_picking") for t in epsilon[0].data))


if __name__ == "__main__":
    unittest.main()
