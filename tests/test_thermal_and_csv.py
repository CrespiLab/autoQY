import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import struct

import numpy as np
import pandas as pd

from autoqy_core.epsilon import (
    EpsilonResult,
    export_epsilon_csv,
    export_epsilon_tsv,
    load_epsilon_table,
)
from autoqy_core.kinetics import _rates
from autoqy_core.io import load_cary_bytes, load_specord_bytes
from autoqy_core.power import load_generic_power_csv
from autoqy_core.smoother import (
    SpectralDataset,
    export_smoothed_text,
    load_spectral_text,
    load_spectral_bytes,
)


class ThermalRateTests(unittest.TestCase):
    def test_both_thermal_directions_are_included(self):
        change = _rates(
            np.array([2.0, 3.0]), 0, np.array([4e-7, 5e-7]), 0, 0,
            np.ones(2), np.ones(2), 1, np.zeros(2), 1,
            0.1, 0.2,
        )
        np.testing.assert_allclose(change, [-0.1, 0.1], atol=1e-15)

    def test_forward_rate_defaults_to_zero_for_compatibility(self):
        change = _rates(
            np.array([2.0, 3.0]), 0, np.array([4e-7, 5e-7]), 0, 0,
            np.ones(2), np.ones(2), 1, np.zeros(2), 1, 0.1,
        )
        np.testing.assert_allclose(change, [0.3, -0.3], atol=1e-15)


class CsvTests(unittest.TestCase):
    def setUp(self):
        values = np.array([1.0, 2.0])
        self.result = EpsilonResult(
            wavelengths=np.array([400.0, 401.0]),
            absorbance=values[:, None],
            concentrations_m=np.array([1e-4]),
            path_lengths_cm=np.array([1.0]),
            individual=(values * 1e4)[:, None],
            mean=values * 1e4,
            standard_deviation=np.array([10.0, 20.0]),
            standard_error=np.array([10.0, 20.0]),
        )

    def test_epsilon_loader_accepts_csv_and_legacy_tsv(self):
        for text in (export_epsilon_csv(self.result, ["sample"]),
                     export_epsilon_tsv(self.result, ["sample"])):
            loaded, labels = load_epsilon_table(text)
            self.assertEqual(labels, ["sample"])
            np.testing.assert_array_equal(loaded.mean, self.result.mean)

    def test_generic_power_csv_uses_milliwatts(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "power.csv"
            pd.DataFrame({"power_mw": [1.2, 1.3]}).to_csv(path, index=False)
            np.testing.assert_array_equal(load_generic_power_csv(path), [1.2, 1.3])

    def test_smoothed_export_defaults_to_csv(self):
        dataset = SpectralDataset(
            np.array([400.0, 401.0]), np.array([0.0]),
            np.array([[0.1], [0.2]]), source_format="spectragryph",
        )
        text = export_smoothed_text(dataset, dataset.absorbance)
        self.assertTrue(text.startswith("Wavelength,0"))

    def test_spectral_csv_accepts_headered_and_headerless_tables(self):
        headered = load_spectral_text(
            "Wavelength,Spectrum\n400,0.1\n401,0.2\n402,0.3\n", "csv"
        )
        headerless = load_spectral_text(
            "400,0.1\n401,0.2\n402,0.3\n", "csv"
        )
        for dataset in (headered, headerless):
            np.testing.assert_array_equal(dataset.wavelengths, [400, 401, 402])
            np.testing.assert_allclose(dataset.absorbance[:, 0], [0.1, 0.2, 0.3])

    def test_returned_csv_text_uses_lf_line_endings(self):
        dataset = SpectralDataset(
            np.array([400.0, 401.0]), np.array([0.0]),
            np.array([[0.1], [0.2]]), source_format="csv",
        )
        for text in (
            export_smoothed_text(dataset, dataset.absorbance),
            export_epsilon_csv(self.result, ["sample"]),
        ):
            self.assertNotIn("\r", text)


class VendorBinaryTests(unittest.TestCase):
    def test_specord_winaspect_binary_loads_multiple_cycles(self):
        wavelengths = np.array([400, 401, 402, 403], dtype="<f4")
        signals = np.array([[0.1, 0.2, 0.3, 0.4],
                            [0.5, 0.6, 0.7, 0.8]], dtype="<f4")
        header = (
            b"[GENERAL]\r\nXUNITS=nm\r\nNPOINTS=4\r\nNCYCL=2\r\n"
            b"[MESS]\r\nORIGIN=SPECORD 200 PLUS\r\n[DATA]\r\nXDATA=\r\n"
        )
        payload = (header + wavelengths.tobytes() + b"\r\nYDATA=\r\n"
                   + signals.tobytes())
        loaded_wavelengths, loaded_signals = load_specord_bytes(payload)
        np.testing.assert_array_equal(loaded_wavelengths, wavelengths)
        np.testing.assert_allclose(loaded_signals, signals.T)
        detected = load_spectral_bytes(payload, "auto")
        self.assertEqual(detected.source_format, "specord")
        self.assertEqual(detected.absorbance.shape, (4, 2))

    def test_cary_winu_v_binary_is_detected_and_ordered(self):
        wavelengths = np.arange(430.0, 399.0, -1.0, dtype=float)
        signals = np.linspace(0.1, 0.4, len(wavelengths))
        header = bytes([17]) + b"Varian UV-VIS-NIR" + bytes(110)
        stream = b"".join(
            struct.pack("<ff", wavelength, signal)
            for wavelength, signal in zip(wavelengths, signals)
        )
        payload = header + stream + bytes(32)
        loaded_wavelengths, loaded_signals = load_cary_bytes(payload)
        np.testing.assert_array_equal(loaded_wavelengths, wavelengths[::-1])
        np.testing.assert_allclose(loaded_signals[:, 0], signals[::-1], rtol=1e-6)
        detected = load_spectral_bytes(payload, "auto")
        self.assertEqual(detected.source_format, "agilent_cary")
        self.assertEqual(detected.absorbance.shape, (31, 1))

    def test_vendor_binary_rejects_wrong_magic(self):
        with self.assertRaisesRegex(ValueError, "Cary"):
            load_cary_bytes(bytes(128))
        with self.assertRaisesRegex(ValueError, "SPECORD"):
            load_specord_bytes(b"not a spectrum")


if __name__ == "__main__":
    unittest.main()
