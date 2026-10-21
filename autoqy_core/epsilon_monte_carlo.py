"""Monte Carlo propagation conditional on fixed PSS reconstruction absorbance.

The model uses truncated Gaussian marginals for reference ε, with a Gaussian
copula estimated from complete replicate curves. Calibrated subtraction reuses
the same sampled reference peak. Physical reconstruction uses a nonnegative
least-squares projection; signed measured PSS absorbance is never changed.
"""

from dataclasses import replace
import warnings

import numpy as np
from scipy.special import ndtr
from scipy.stats import binom, truncnorm

from .epsilon_uncertainty import (
    EpsilonUncertaintySummary, _absorbance_residual_rmse, _column,
    _read_epsilon_frame, load_epsilon_envelope,
)


SCOPE = ("Reference ε, NMR composition, irradiation input, and kinetic-fit residual "
         "contributions; PSS reconstruction absorbance is fixed. No separate PSS ε is required.")


def _positive_normal(mean, sigma, probability):
    """Draw conditional nonnegative values, without clipping Gaussian draws to zero."""
    mean, sigma, probability = np.broadcast_arrays(mean, sigma, probability)
    if np.any(~np.isfinite(mean)) or np.any(~np.isfinite(sigma)) or np.any(sigma < 0):
        raise ValueError("Finite reference means and nonnegative uncertainties are required")
    if np.any((sigma == 0) & (mean < 0)):
        raise ValueError("A negative reference ε has zero uncertainty; review the reference export")
    result = mean.copy().astype(float)
    active = sigma > 0
    result[active] = truncnorm.ppf(
        np.clip(probability[active], np.finfo(float).eps, 1 - np.finfo(float).eps),
        -mean[active] / sigma[active], np.inf, loc=mean[active], scale=sigma[active],
    )
    return result


def _reference_sampler(frame, metric):
    mean = _column(frame, "Epsilon_mean_M-1_cm-1", "reference")
    sigma = _column(frame, "Epsilon_SD_M-1_cm-1" if metric == "sd"
                    else "Epsilon_SEM_M-1_cm-1", "reference")
    columns = [name for name in frame if name.startswith("Epsilon_M-1_cm-1__")]
    if len(columns) < 2 and np.any(sigma > 0):
        raise ValueError("Monte Carlo requires at least two complete reference ε replicate curves. "
                         "Use ε ranges for older exports without replicates.")
    if columns:
        replicates = frame[columns].to_numpy(float)
        if not np.isfinite(replicates).all():
            raise ValueError("Reference ε replicate curves must be finite")
        deviations = replicates - replicates.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(deviations, axis=1)
        factors = np.divide(deviations, norm[:, None], out=np.zeros_like(deviations),
                            where=norm[:, None] > 0)
        if np.any((norm == 0) & (sigma > 0)):
            raise ValueError("Reference errors disagree with its identical replicate curves")
    else:
        factors = np.zeros((len(mean), 1))

    def sample(rng):
        z = factors @ rng.normal(size=factors.shape[1])
        return _positive_normal(mean, sigma, ndtr(z))
    # Validate deterministic contradictions before fitting any simulations.
    _positive_normal(mean, sigma, np.full(mean.shape, .5))
    return sample


def _scalar(frame, name):
    values = _column(frame, name, "NMR reconstruction")
    if not np.allclose(values, values[0], rtol=1e-7, atol=1e-10):
        raise ValueError(f"{name} must be constant across the reconstruction")
    return float(values[0])


def _reconstruction_sampler(frame, reference_frame, reference_axis):
    required = {"Reactant_normalized", "PSS_normalized_to_reactant",
                "Reactant_reconstructed_M-1_cm-1", "Product_epsilon_raw_M-1_cm-1",
                "Composition_uncertainty_percentage_points"}
    if not required.issubset(frame):
        raise ValueError("Monte Carlo needs the calibrated reconstruction and NMR composition "
                         "uncertainty. Re-export it from Spectral Treatment, or choose ε ranges.")
    x = _column(frame, "Wavelength_nm", "NMR reconstruction")
    if x[0] < reference_axis[0] or x[-1] > reference_axis[-1]:
        raise ValueError("The shared reference does not cover the reconstructed spectrum")
    stable = _column(frame, "Reactant_normalized", "NMR reconstruction")
    mixture = _column(frame, "PSS_normalized_to_reactant", "NMR reconstruction")
    calibrated = _column(frame, "Reactant_reconstructed_M-1_cm-1", "NMR reconstruction")
    raw = _column(frame, "Product_epsilon_raw_M-1_cm-1", "NMR reconstruction")
    mask = np.abs(stable) > .1 * np.max(np.abs(stable))
    scale = float(np.median(calibrated[mask] / stable[mask]))
    reference = np.interp(x, reference_axis,
                          reference_frame["Epsilon_mean_M-1_cm-1"].to_numpy(float))
    if "Calibration_reference_epsilon_M-1_cm-1" in frame:
        saved_reference = _column(frame, "Calibration_reference_epsilon_M-1_cm-1", "NMR reconstruction")
        if not np.allclose(saved_reference, reference, rtol=2e-5, atol=1e-5):
            raise ValueError("The supplied reference differs from the shared reference saved in the reconstruction")
    if scale <= 0 or not np.isclose(scale, np.max(reference), rtol=2e-5, atol=1e-6):
        raise ValueError("The reconstructed spectrum cannot be matched to the supplied shared "
                         "reference calibration. Choose ε ranges or re-export the matching pair.")
    if not np.allclose(calibrated, stable * scale, rtol=2e-5, atol=1e-5):
        raise ValueError("The stored calibration does not use a shared reference scale")
    a, b = raw / scale - mixture, raw / scale - stable
    if "Stable_fraction_in_PSS" in frame:
        fraction = _scalar(frame, "Stable_fraction_in_PSS")
    else:
        denominator = float(b @ b)
        if denominator <= np.finfo(float).eps:
            raise ValueError("Cannot recover the NMR composition from this older reconstruction")
        fraction = float(a @ b / denominator)
    if not 0 <= fraction < 1 or not np.allclose(
            (mixture - fraction * stable) * scale / (1 - fraction), raw,
            rtol=2e-5, atol=2e-4):
        raise ValueError("Cannot validate the NMR composition in this reconstruction")
    sigma = _scalar(frame, "Composition_uncertainty_percentage_points") / 100
    if sigma < 0:
        raise ValueError("NMR composition uncertainty must be nonnegative")

    def sample(rng, sampled_reference):
        f = (fraction if sigma == 0 else float(truncnorm.rvs(
            -fraction / sigma, (1 - fraction) / sigma,
            loc=fraction, scale=sigma, random_state=rng)))
        sampled_scale = np.max(np.interp(x, reference_axis, sampled_reference))
        raw_curve = sampled_scale * (mixture - f * stable) / (1 - f)
        # Exact solution of nonnegative least-squares subtraction at fixed f.
        # This is the constrained reconstruction estimator, not clipped ε noise.
        return np.maximum(raw_curve, 0)
    return sample, fraction


def paired_epsilon_sampler(reactant, product, metric="sd"):
    """Build the joint sampler before starting costly kinetic fits."""
    frames = [_read_epsilon_frame(e.source_path) for e in (reactant, product)]
    derived = [i for i, f in enumerate(frames) if "Product_epsilon_M-1_cm-1" in f]
    assumptions = [
        "Reference ε: nonnegative truncated Gaussian marginals; Gaussian copula correlations "
        "estimated from whole replicate deviations. SD/SEM describe the untruncated measurement model.",
        "Kinetic noise: centred whole-spectrum residual bootstrap across timestamps; "
        "residual spectra assumed exchangeable. Local optimizer errors are not added again.",
        "Power/photon-flux and NMR error fields are interpreted as one standard uncertainty; "
        "their errors, reference repeatability, and kinetic residuals are assumed independent.",
    ]
    if not derived:
        samplers = [_reference_sampler(f, metric) for f in frames]
        assumptions.append("Separately supplied reference spectra are assumed independent.")
        return lambda rng: tuple(s(rng) for s in samplers), tuple(assumptions)
    if len(derived) != 1:
        raise ValueError("Monte Carlo needs one measured reference and one matching calibrated "
                         "reconstruction, or two measured references. Choose ε ranges for this pair.")
    index = derived[0]
    reference_index = 1 - index
    reference_envelope = (reactant, product)[reference_index]
    if "Calibration_reference_epsilon_M-1_cm-1" not in frames[index]:
        assumptions.append("Legacy shared-reference association inferred from the saved calibration "
                           "scale and supplied reference; the original reference identity is not stored.")
    reference_sample = _reference_sampler(frames[reference_index], metric)
    reconstruction_sample, fraction = _reconstruction_sampler(
        frames[index], frames[reference_index], reference_envelope.wavelengths)
    assumptions.append(f"Shared calibration reused in both species; stable fraction in "
                       f"reconstruction PSS = {fraction:.10g}. Fixed normalized absorbance "
                       "curves; nonnegative least-squares subtraction for each draw.")

    def sample(rng):
        reference = reference_sample(rng)
        reconstructed = reconstruction_sample(rng, reference)
        return ((reconstructed, reference) if index == 0 else (reference, reconstructed))
    return sample, tuple(assumptions)


def _endpoint_error(samples, coverage, yield_bounds=(0, 1)):
    """Conservative endpoint sampling error from batches and binomial order statistics.

    Small batches can miss the tails and underestimate quantile error. Use the
    larger of that estimate and a distribution-free 95% quantile bracket scaled
    to an approximate standard error. Known yield bounds cover unresolved tails.
    """
    samples = np.asarray(samples)
    tail = (1 - coverage) / 2
    groups = np.array_split(samples, 4)
    intervals = np.asarray([np.quantile(g, [tail, 1 - tail], axis=0) for g in groups])
    batch_error = np.std(intervals, axis=0, ddof=1) / 2
    ordered = np.sort(samples, axis=0)
    errors = []
    for probability in (tail, 1 - tail):
        low_rank = int(binom.ppf(.025, len(samples), probability))
        high_rank = int(binom.ppf(.975, len(samples), probability)) + 1
        lower = (np.full(samples.shape[1], yield_bounds[0]) if low_rank == 0
                 else ordered[low_rank - 1])
        upper = (np.full(samples.shape[1], yield_bounds[1]) if high_rank > len(samples)
                 else ordered[high_rank - 1])
        errors.append((upper - lower) / (2 * 1.96))
    result = np.maximum(batch_error, errors)
    result[:, np.ptp(samples, axis=0) == 0] = 0
    return result


def run_monte_carlo(data, reactant, product, *, error_metric="sd", samples=256,
                    seed=0, coverage=.95, tolerance=.005):
    from .pipeline import run_analysis_pipeline

    sample_pair, assumptions = paired_epsilon_sampler(reactant, product, error_metric)
    rng = np.random.default_rng(seed)
    nominal = run_analysis_pipeline(data)
    if not nominal.yield_fit.optimizer_success:
        raise ValueError("The nominal kinetic fit failed; Monte Carlo cannot start")
    model = nominal.yield_fit.concentrations @ np.vstack((nominal.epsilon_r, nominal.epsilon_p))
    model *= data.path_length_cm
    if nominal.yield_fit.absorbance_correction is not None:
        model += nominal.yield_fit.absorbance_correction
    residuals = nominal.absorbance.T - model
    residuals -= residuals.mean(axis=0, keepdims=True)
    rows = np.searchsorted(data.wavelengths, nominal.wavelengths)
    records, failures = [], []
    attempted = 0
    endpoint_error = np.full((2, 2), np.inf)
    for iteration in range(samples):
        attempted += 1
        try:
            er, ep = sample_pair(rng)
            absorbance = data.absorbance.copy()
            bootstrap = residuals[rng.integers(0, len(residuals), size=len(residuals))]
            absorbance[rows] = (model + bootstrap).T
            power, flux = data.power_mw, data.photon_flux_mol_s
            if flux is None:
                power = float(_positive_normal(power, data.power_error_mw, rng.random()))
            else:
                flux = float(_positive_normal(flux, data.photon_flux_error_mol_s, rng.random()))
            trial_data = replace(
                data, epsilon_r=(reactant.wavelengths, er), epsilon_p=(product.wavelengths, ep),
                absorbance=absorbance, power_mw=power, photon_flux_mol_s=flux,
                power_error_mw=0, photon_flux_error_mol_s=0,
            )
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Initial product is .*", category=RuntimeWarning)
                trial = run_analysis_pipeline(trial_data, include_power_variations=False)
            fit = trial.yield_fit
            if not fit.optimizer_success or not np.isfinite(fit.values).all():
                raise ValueError("Kinetic optimizer failed or returned non-finite yields")
            total = fit.concentrations.sum(axis=1)
            fraction = np.divide(fit.concentrations[:, 0], total,
                                 out=np.zeros_like(total), where=total != 0)
            residual = trial.concentration_fit.fractions[:, 0] - fraction
            records.append((fit.values, fit.standard_errors, fit.active_bounds,
                            fit.jacobian_condition, trial.concentration_fit.concentrations,
                            fit.concentrations, residual,
                            _absorbance_residual_rmse(trial, data.path_length_cm)))
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as error:
            failures.append(f"Simulation {iteration + 1}: {error}")
        if len(records) >= 128 and attempted % 64 == 0:
            endpoint_error = _endpoint_error([r[0] for r in records], coverage, data.yield_bounds)
            if not failures and np.all(2 * endpoint_error <= tolerance):
                break
    if len(records) < max(16, attempted // 2):
        raise ValueError(f"Only {len(records)} of {attempted} Monte Carlo fits succeeded. "
                         + (failures[0] if failures else "Increase the number of simulations."))
    yields = np.asarray([r[0] for r in records])
    if len(records) >= 16:
        endpoint_error = _endpoint_error(yields, coverage, data.yield_bounds)
    tail = (1 - coverage) / 2
    lower, upper = np.quantile(yields, [tail, 1 - tail], axis=0)
    std = np.std(yields, axis=0, ddof=1)
    bounds = np.asarray([r[2] for r in records], bool)
    cd = np.asarray([r[4] for r in records])
    cf = np.asarray([r[5] for r in records])
    fr = np.asarray([r[6] for r in records])
    rmse = np.asarray([r[7] for r in records])
    summary = EpsilonUncertaintySummary(
        method="monte_carlo", error_metric=error_metric,
        bound_combination_count=len(records),
        reactant_source_schema=reactant.source_schema, product_source_schema=product.source_schema,
        reactant_source_path=reactant.source_path, product_source_path=product.source_path,
        reactant_error_metric=reactant.error_metric, product_error_metric=product.error_metric,
        constrained_negative_points=(reactant.constrained_negative_points, product.constrained_negative_points),
        optimizer_power_errors=nominal.yield_errors, epsilon_errors=np.zeros(2), combined_errors=std,
        epsilon_yield_minimum=lower, epsilon_yield_maximum=upper,
        concentration_data_minimum=np.quantile(cd, tail, axis=0),
        concentration_data_maximum=np.quantile(cd, 1 - tail, axis=0),
        concentration_fit_minimum=np.quantile(cf, tail, axis=0),
        concentration_fit_maximum=np.quantile(cf, 1 - tail, axis=0),
        fraction_residual_minimum=np.quantile(fr, tail, axis=0),
        fraction_residual_maximum=np.quantile(fr, 1 - tail, axis=0),
        absorbance_residual_rmse_minimum=float(rmse.min()),
        absorbance_residual_rmse_maximum=float(rmse.max()),
        bound_labels=tuple(f"Simulation {i + 1}" for i in range(len(records))),
        bound_yields=yields, bound_optimizer_errors=np.asarray([r[1] for r in records]),
        bound_fraction_rmse=np.sqrt(np.mean(fr ** 2, axis=1)), bound_absorbance_rmse=rmse,
        bound_active_bounds=bounds, bound_jacobian_conditions=np.asarray([r[3] for r in records]),
        bound_optimizer_success=tuple(True for _ in records),
        coverage_probability=coverage, interval_lower=lower, interval_upper=upper,
        sample_mean=yields.mean(axis=0), sample_standard_deviation=std,
        attempted_samples=attempted, requested_samples=samples, failed_samples=len(failures),
        boundary_fractions=bounds.mean(axis=0), random_seed=seed,
        converged=bool(not failures and len(records) >= 128 and np.all(2 * endpoint_error <= tolerance)),
        endpoint_mc_error=endpoint_error, scope=SCOPE, sampling_assumptions=assumptions,
        failure_messages=tuple(failures),
    )
    return replace(nominal, yield_errors=std, epsilon_uncertainty=summary)
