"""Create headless AutoQY result figures."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from .output import format_value_uncertainty


def _nipe_window_analysis(result):
    nipe = getattr(result.yield_fit, "nipe", None)
    return getattr(nipe, "window_analysis", None)


def _nipe_fit_display_window(result, timestamps):
    times = np.asarray(timestamps, dtype=float)
    mask = np.ones(times.shape, dtype=bool)
    if result.fit_method != "nipe" or times.size == 0:
        return mask, None
    window = _nipe_window_analysis(result)
    end_relative = getattr(window, "analysis_end_time_s", None)
    if end_relative is None or not np.isfinite(end_relative):
        return mask, None
    end_absolute = min(float(times[-1]), float(times[0]) + float(end_relative))
    tolerance = max(1.0, abs(end_absolute)) * np.finfo(float).eps * 16
    mask = times <= end_absolute + tolerance
    return (mask, float(times[mask][-1])) if np.any(mask) else (
        np.ones(times.shape, dtype=bool), None
    )


def _effective_ab_model_level(result):
    uncertainty = getattr(result, "epsilon_uncertainty", None)
    if uncertainty is not None:
        return getattr(uncertainty, "ab_model_level", None)
    diagnostic = getattr(result, "ab_model_diagnostic", None)
    return getattr(diagnostic, "level", None)


def _quantum_yield_annotation(result):
    phi_rp = "$\\Phi_{\\mathrm{R}\\rightarrow\\mathrm{P}}$"
    phi_pr = "$\\Phi_{\\mathrm{P}\\rightarrow\\mathrm{R}}$"

    def formatted_lines(label, values, errors):
        formatted = [
            format_value_uncertainty(value, error, two_digit_threshold=2)
            for value, error in zip(np.asarray(values) * 100, np.asarray(errors) * 100)
        ]
        values_text = (
            f"{phi_rp}: {formatted[0][0]} ± {formatted[0][1]}%\n"
            f"{phi_pr}: {formatted[1][0]} ± {formatted[1][1]}%"
        )
        return f"{label}\n{values_text}" if label else values_text

    complete = formatted_lines(
        None, result.yield_fit.values, result.yield_errors,
    )
    if result.fit_method != "nipe" or _effective_ab_model_level(result) != "stop":
        return complete

    window = _nipe_window_analysis(result)
    if window is None:
        return complete + "\n\nRed NIPE flag: no reliable pre-plateau estimate available"
    errors = window.extrapolated_standard_errors
    uncertainty = getattr(result, "epsilon_uncertainty", None)
    includes_epsilon = (
        uncertainty is not None
        and getattr(uncertainty, "nipe_window_combined_errors", None) is not None
    )
    if includes_epsilon:
        errors = uncertainty.nipe_window_combined_errors
    label = "Recommended pre-plateau NIPE estimate"
    if includes_epsilon:
        label += " (includes ε range)"
    return (
        formatted_lines(
            "Complete-trace NIPE (red flag)",
            result.yield_fit.values,
            result.yield_errors,
        )
        + "\n\n"
        + formatted_lines(label, window.extrapolated_values, errors)
    )


def write_figure(path, result, data, residual_percentile=100):
    if not 0 < residual_percentile <= 100:
        raise ValueError("residual_percentile must be greater than 0 and at most 100")
    times = np.asarray(data.timestamps)
    measured = result.concentration_fit.concentrations
    fitted = result.yield_fit.concentrations
    fitted_fraction = fitted[:, 0] / fitted.sum(axis=1)
    fraction_residual = result.concentration_fit.fractions[:, 0] - fitted_fraction
    epsilon = np.vstack((result.epsilon_r, result.epsilon_p))
    measured_absorbance = result.absorbance.T
    fitted_absorbance = fitted @ epsilon * data.path_length_cm
    if result.yield_fit.absorbance_correction is not None:
        fitted_absorbance = fitted_absorbance + result.yield_fit.absorbance_correction
    absorbance_residual = measured_absorbance - fitted_absorbance
    uncertainty = result.epsilon_uncertainty
    fit_mask, fit_end_time = _nipe_fit_display_window(result, times)
    fit_times = times[fit_mask]

    blue, orange, brown = "#346aa9", "#e16203", "#8a6642"
    figure, axes = plt.subplots(
        2, 2, figsize=(13, 7), constrained_layout=True,
        gridspec_kw={"wspace": 0.16},
    )
    concentration, spectra, residual, heatmap = axes.flat

    if uncertainty is not None:
        for index, colour in enumerate((blue, orange)):
            concentration.errorbar(
                times, measured[:, index],
                yerr=np.vstack((
                    measured[:, index] - uncertainty.concentration_data_minimum[:, index],
                    uncertainty.concentration_data_maximum[:, index] - measured[:, index],
                )),
                fmt="none", ecolor=colour, elinewidth=0.8, alpha=0.35, capsize=2,
            )
            concentration.fill_between(
                fit_times, uncertainty.concentration_fit_minimum[fit_mask, index],
                uncertainty.concentration_fit_maximum[fit_mask, index],
                color=colour, alpha=0.13,
            )
    nominal = " (nominal ε)" if uncertainty is not None else ""
    concentration.scatter(times, measured[:, 0], s=24, facecolors="none",
                          edgecolors=blue, label=f"Reactant data{nominal}")
    concentration.scatter(times, measured[:, 1], s=24, facecolors="none",
                          edgecolors=orange, label=f"Product data{nominal}")
    concentration.plot(fit_times, fitted[fit_mask, 0], color=blue, linewidth=2,
                       zorder=4, label=f"Reactant fit{nominal}")
    concentration.plot(fit_times, fitted[fit_mask, 1], color=orange, linewidth=2,
                       zorder=4, label=f"Product fit{nominal}")
    if fit_end_time is not None:
        concentration.axvline(fit_end_time, color=brown, linestyle="--", linewidth=1.4)
        concentration.annotate(
            "NIPE fit window ends", xy=(fit_end_time, 1),
            xycoords=("data", "axes fraction"), xytext=(4, -4),
            textcoords="offset points", ha="left", va="top", fontsize=8, color=brown,
        )
    concentration.set(title=(
                          "NIPE concentrations: fit limited to accepted window"
                          if result.fit_method == "nipe" else
                          "Concentrations: nominal ε and ε-bound ranges"
                          if uncertainty is not None else "Concentration fit"),
                      xlabel="Irradiation time (s)",
                      ylabel="Concentration (mol/L)")
    concentration.legend(frameon=False, ncol=2)

    if uncertainty is not None:
        residual.fill_between(
            times, uncertainty.fraction_residual_minimum,
            uncertainty.fraction_residual_maximum,
            color=blue, alpha=0.16, label="ε-bound range",
        )
    residual.plot(times, fraction_residual, "o-", color=blue, markersize=4,
                  label=f"Nominal ε residual" if uncertainty is not None else None)
    residual.axhline(0, color="black", linewidth=0.8)
    residual.set(title="Reactant fraction residual", xlabel="Irradiation time (s)",
                 ylabel="Fraction data - fit")
    if uncertainty is not None:
        residual.legend(frameon=False)

    normalization = Normalize(times.min(), times.max())
    colour_map = plt.get_cmap("RdBu_r")
    for time, spectrum in zip(times, measured_absorbance):
        spectra.plot(result.wavelengths, spectrum, color=colour_map(normalization(time)),
                     linewidth=1)
    spectra.text(0.98, 0.96,
                 _quantum_yield_annotation(result),
                 transform=spectra.transAxes, ha="right", va="top", fontsize=8.5,
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
    heatmap_title = "Absorbance residuals"
    if uncertainty is not None:
        heatmap_title = (
            "Absorbance residuals (nominal ε)\n"
            f"ε-bound RMSE {uncertainty.absorbance_residual_rmse_minimum:.3g}–"
            f"{uncertainty.absorbance_residual_rmse_maximum:.3g}"
        )
    heatmap.set(title=heatmap_title, xlabel="Wavelength (nm)",
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
