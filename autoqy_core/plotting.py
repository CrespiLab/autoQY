"""Create headless AutoQY result figures."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from .output import format_value_uncertainty


def write_figure(path, result, data, residual_percentile=100):
    if not 0 < residual_percentile <= 100:
        raise ValueError("residual_percentile must be greater than 0 and at most 100")
    times = data.timestamps
    measured = result.concentration_fit.concentrations
    fitted = result.yield_fit.concentrations
    fitted_fraction = fitted[:, 0] / fitted.sum(axis=1)
    fraction_residual = result.concentration_fit.fractions[:, 0] - fitted_fraction
    epsilon = np.vstack((result.epsilon_r, result.epsilon_p))
    measured_absorbance = result.absorbance.T
    fitted_absorbance = fitted @ epsilon * data.path_length_cm
    absorbance_residual = measured_absorbance - fitted_absorbance

    blue, orange = "#346aa9", "#e16203"
    figure, axes = plt.subplots(
        2, 2, figsize=(13, 7), constrained_layout=True,
        gridspec_kw={"wspace": 0.16},
    )
    concentration, spectra, residual, heatmap = axes.flat

    concentration.scatter(times, measured[:, 0], s=24, facecolors="none",
                          edgecolors=blue, label="Reactant data")
    concentration.scatter(times, measured[:, 1], s=24, facecolors="none",
                          edgecolors=orange, label="Product data")
    concentration.plot(times, fitted[:, 0], color=blue, label="Reactant fit")
    concentration.plot(times, fitted[:, 1], color=orange, label="Product fit")
    concentration.set(title="Concentration fit", xlabel="Irradiation time (s)",
                      ylabel="Concentration (mol/L)")
    concentration.legend(frameon=False, ncol=2)

    residual.plot(times, fraction_residual, "o-", color=blue, markersize=4)
    residual.axhline(0, color="black", linewidth=0.8)
    residual.set(title="Reactant fraction residual", xlabel="Irradiation time (s)",
                 ylabel="Fraction data - fit")

    normalization = Normalize(times.min(), times.max())
    colour_map = plt.get_cmap("RdBu_r")
    for time, spectrum in zip(times, measured_absorbance):
        spectra.plot(result.wavelengths, spectrum, color=colour_map(normalization(time)),
                     linewidth=1)
    formatted = [format_value_uncertainty(value, error) for value, error in
                 zip(result.yield_fit.values * 100, result.yield_errors * 100)]
    phi_rp = "$\\Phi_{\\mathrm{R}\\rightarrow\\mathrm{P}}$"
    phi_pr = "$\\Phi_{\\mathrm{P}\\rightarrow\\mathrm{R}}$"
    spectra.text(0.98, 0.96,
                 f"{phi_rp}: {formatted[0][0]} \u00b1 {formatted[0][1]}%\n"
                 f"{phi_pr}: {formatted[1][0]} \u00b1 {formatted[1][1]}%",
                 transform=spectra.transAxes, ha="right", va="top",
                 bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
    spectra.set(title="Absorption spectra over time", xlabel="Wavelength (nm)",
                ylabel="Absorbance")
    wavelength_limits = result.wavelengths[0], result.wavelengths[-1]
    spectra.set_xlim(*wavelength_limits)
    spectra.margins(x=0)
    figure.colorbar(plt.cm.ScalarMappable(normalization, colour_map), ax=spectra,
                    label="Irradiation time (s)")

    limit = np.percentile(np.abs(absorbance_residual), residual_percentile)
    limit = limit or np.finfo(float).eps
    image = heatmap.imshow(absorbance_residual, aspect="auto", cmap=colour_map,
                           vmin=-limit, vmax=limit, origin="upper",
                           extent=(*wavelength_limits, times[-1], times[0]))
    heatmap.set(title="Absorbance residuals", xlabel="Wavelength (nm)",
                ylabel="Irradiation time (s)")
    heatmap.set_xlim(*wavelength_limits)
    figure.colorbar(image, ax=heatmap, label="Data - fit")

    for axis in (concentration, residual, spectra):
        axis.grid(axis="y", alpha=0.2)
    path = Path(path)
    stem = path.with_suffix("") if path.suffix.lower() in {".png", ".svg"} else path
    png_path, svg_path = stem.with_suffix(".png"), stem.with_suffix(".svg")
    figure.savefig(png_path, dpi=300, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, svg_path
