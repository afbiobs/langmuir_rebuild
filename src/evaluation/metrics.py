"""
Evaluation metrics for Langmuir circulation spacing predictions.

All functions operate on arrays of predicted and observed spacings in metres.
Every function documents units explicitly.

Anti-saturation check: attractor_test and saturation_audit are the primary
structural health checks. They must be run after every module is completed.
"""

import numpy as np
import pandas as pd
from typing import Callable


def _validate_pairwise_inputs(
    predicted: np.ndarray,
    observed: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate and coerce paired predicted/observed spacing arrays [m]."""
    predicted = np.asarray(predicted, dtype=float)
    observed = np.asarray(observed, dtype=float)
    if len(predicted) != len(observed):
        raise ValueError(
            f"Length mismatch: predicted={len(predicted)}, observed={len(observed)}"
        )
    return predicted, observed


def spacing_rmse(predicted: np.ndarray, observed: np.ndarray) -> float:
    """
    Root-mean-square error of spacing predictions.

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]

    Returns:
        RMSE [m]
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    return float(np.sqrt(np.mean((predicted - observed) ** 2)))


def spacing_mae(predicted: np.ndarray, observed: np.ndarray) -> float:
    """
    Mean absolute error of spacing predictions.

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]

    Returns:
        MAE [m]
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    return float(np.mean(np.abs(predicted - observed)))


def spacing_rmse_ratio(predicted: np.ndarray, observed: np.ndarray) -> float:
    """
    Root-mean-square relative spacing error [%].

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]

    Returns:
        RMSE of (predicted - observed) / observed [%]
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    if np.any(observed <= 0.0):
        raise ValueError("observed spacings must be positive for relative metrics.")
    residual_ratio = (predicted - observed) / observed
    return float(100.0 * np.sqrt(np.mean(residual_ratio ** 2)))


def spacing_mae_ratio(predicted: np.ndarray, observed: np.ndarray) -> float:
    """
    Mean absolute relative spacing error [%].

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]

    Returns:
        MAE of |predicted - observed| / observed [%]
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    if np.any(observed <= 0.0):
        raise ValueError("observed spacings must be positive for relative metrics.")
    residual_ratio = np.abs(predicted - observed) / observed
    return float(100.0 * np.mean(residual_ratio))


def spacing_bias(predicted: np.ndarray, observed: np.ndarray) -> float:
    """
    Mean signed spacing error [m].

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]

    Returns:
        Mean(predicted - observed) [m]
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    return float(np.mean(predicted - observed))


def pearson_correlation(predicted: np.ndarray, observed: np.ndarray) -> float:
    """
    Pearson linear correlation coefficient between prediction and observation [-].

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]

    Returns:
        Pearson r [-], or nan when one series is constant
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    if len(predicted) < 2:
        return float("nan")
    if np.std(predicted) == 0.0 or np.std(observed) == 0.0:
        return float("nan")
    return float(np.corrcoef(predicted, observed)[0, 1])


def spearman_correlation(predicted: np.ndarray, observed: np.ndarray) -> float:
    """
    Spearman rank correlation coefficient between prediction and observation [-].

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]

    Returns:
        Spearman ρ [-], or nan when one ranked series is constant
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    if len(predicted) < 2:
        return float("nan")
    ranked_predicted = pd.Series(predicted).rank(method="average").to_numpy(dtype=float)
    ranked_observed = pd.Series(observed).rank(method="average").to_numpy(dtype=float)
    if np.std(ranked_predicted) == 0.0 or np.std(ranked_observed) == 0.0:
        return float("nan")
    return float(np.corrcoef(ranked_predicted, ranked_observed)[0, 1])


def hit_rate_within_fraction(
    predicted: np.ndarray,
    observed: np.ndarray,
    tolerance: float = 0.30,
) -> float:
    """
    Fraction of cases predicted within ±tolerance of the observed spacing [-].

    Parameters:
        predicted: Predicted windrow spacings [m]
        observed:  Observed windrow spacings [m]
        tolerance: Fractional hit tolerance [-]

    Returns:
        Hit rate [-]
    """
    predicted, observed = _validate_pairwise_inputs(predicted, observed)
    if np.any(observed <= 0.0):
        raise ValueError("observed spacings must be positive for hit-rate metrics.")
    return float(np.mean(np.abs(predicted - observed) <= tolerance * observed))


def range_coverage(
    predicted: np.ndarray,
    observed_range: tuple[float, float],
) -> float:
    """
    Fraction of the observed range spanned by the predictions [-].

    Parameters:
        predicted:       Predicted windrow spacings [m]
        observed_range:  (min, max) observed spacing [m]

    Returns:
        (max(predicted) - min(predicted)) / (obs_max - obs_min) [-]
    """
    predicted = np.asarray(predicted, dtype=float)
    obs_min, obs_max = float(observed_range[0]), float(observed_range[1])
    obs_span = obs_max - obs_min
    if obs_span <= 0.0:
        return float("nan")
    return float((np.max(predicted) - np.min(predicted)) / obs_span)


def tail_coverage(
    predicted: np.ndarray,
    observed: np.ndarray,
    low_threshold: float = 60.0,
    high_threshold: float = 120.0,
    tolerance: float = 0.30,
) -> dict:
    """
    Fraction of tail observations predicted within ±tolerance of observed.

    Tail observations are those where observed spacing is outside the central
    range [low_threshold, high_threshold]. These are the physically interesting
    cases: near-onset cells (narrow spacing) and coarsened cells (wide spacing).

    Parameters:
        predicted:       Predicted spacings [m]
        observed:        Observed spacings [m]
        low_threshold:   Upper bound of low-spacing tail [m] (default: 60)
        high_threshold:  Lower bound of high-spacing tail [m] (default: 120)
        tolerance:       Fractional tolerance for "hit" (default: 0.30 = ±30%)

    Returns:
        dict with keys:
            low_tail_coverage:  fraction of low-tail obs predicted within ±tolerance [-]
            high_tail_coverage: fraction of high-tail obs predicted within ±tolerance [-]
            overall_tail_coverage: combined fraction [-]
            n_low_tail:   number of low-tail observations [-]
            n_high_tail:  number of high-tail observations [-]
    """
    predicted = np.asarray(predicted, dtype=float)
    observed = np.asarray(observed, dtype=float)

    low_mask = observed < low_threshold
    high_mask = observed > high_threshold

    def _coverage(pred, obs):
        if len(obs) == 0:
            return float("nan"), 0
        hits = np.abs(pred - obs) <= tolerance * obs
        return float(np.mean(hits)), len(obs)

    low_cov, n_low = _coverage(predicted[low_mask], observed[low_mask])
    high_cov, n_high = _coverage(predicted[high_mask], observed[high_mask])

    n_tail = n_low + n_high
    if n_tail == 0:
        overall = float("nan")
    else:
        n_low_hits = int(round(low_cov * n_low)) if not np.isnan(low_cov) else 0
        n_high_hits = int(round(high_cov * n_high)) if not np.isnan(high_cov) else 0
        overall = (n_low_hits + n_high_hits) / n_tail

    return {
        "low_tail_coverage": low_cov,
        "high_tail_coverage": high_cov,
        "overall_tail_coverage": overall,
        "n_low_tail": n_low,
        "n_high_tail": n_high,
    }


def dynamic_range(predicted: np.ndarray) -> float:
    """
    Ratio of 90th to 10th percentile of predictions.

    A value near 1.0 indicates the predictions are clustered (attractor behaviour).
    The observed spacing distribution has a dynamic range of ~518/37 ≈ 14 over the
    full dataset; the 90th/10th ratio for the real data is approximately 3–5.

    Parameters:
        predicted: Predicted spacings [m]

    Returns:
        p90 / p10 ratio [-]. Returns nan if p10 == 0.
    """
    predicted = np.asarray(predicted, dtype=float)
    p10 = float(np.percentile(predicted, 10))
    p90 = float(np.percentile(predicted, 90))
    if p10 == 0.0:
        return float("nan")
    return p90 / p10


def attractor_test(
    predicted: np.ndarray,
    observed_range: tuple,
    band_fraction: float = 0.20,
    concentration_threshold: float = 0.50,
) -> dict:
    """
    FAIL if more than concentration_threshold of predictions fall in a band spanning
    less than band_fraction of the observed range.

    This is the primary structural health check for the attractor failure mode:
    a model that predicts ~80 m for all inputs regardless of forcing. Such a model
    can achieve low RMSE on a dataset centred on 80 m while being physically useless.

    Parameters:
        predicted:               Predicted spacings [m]
        observed_range:          (min, max) of observed spacings [m]
        band_fraction:           Width of detection band as fraction of observed range [-]
                                 Default: 0.20 (flag if >50% in <20% of range)
        concentration_threshold: Fraction of predictions to trigger failure [-]
                                 Default: 0.50

    Returns:
        dict with keys:
            passed:        bool — True if attractor NOT detected
            max_fraction:  highest fraction found in any band of width band_fraction [-]
            band_width:    absolute width of detection band [m]
            message:       human-readable diagnostic
    """
    predicted = np.asarray(predicted, dtype=float)
    obs_min, obs_max = float(observed_range[0]), float(observed_range[1])
    obs_span = obs_max - obs_min
    band_width = band_fraction * obs_span

    # Slide a window of width band_width across the observed range
    n_steps = 200
    centres = np.linspace(obs_min + band_width / 2, obs_max - band_width / 2, n_steps)
    max_frac = 0.0
    for c in centres:
        lo, hi = c - band_width / 2, c + band_width / 2
        frac = float(np.mean((predicted >= lo) & (predicted <= hi)))
        if frac > max_frac:
            max_frac = frac

    passed = max_frac <= concentration_threshold
    msg = (
        f"OK: max fraction in any {band_fraction*100:.0f}% band = {max_frac:.2%}"
        if passed
        else (
            f"FAIL: {max_frac:.2%} of predictions in a band spanning "
            f"{band_fraction*100:.0f}% of observed range [{obs_min:.0f}, {obs_max:.0f}] m. "
            f"Attractor behaviour detected."
        )
    )

    return {
        "passed": passed,
        "max_fraction": max_frac,
        "band_width": band_width,
        "message": msg,
    }


def saturation_audit(
    func: Callable,
    input_ranges: dict,
    n_points: int = 50,
    threshold_variation: float = 0.05,
    threshold_fraction: float = 0.50,
) -> dict:
    """
    Sweep func over input ranges. Flag any scalar output that varies less than
    threshold_variation (5%) over more than threshold_fraction (50%) of the input range.

    This implements the saturation rule from AGENTS.md §2. A function whose output
    varies by <5% over >50% of its physically plausible input range is saturated and
    must be fixed before use.

    Parameters:
        func:                Callable accepting keyword arguments from input_ranges.
                             Must return a scalar or dict of scalars.
        input_ranges:        dict mapping parameter name → (min, max) tuple.
                             Only one parameter is varied at a time; others are held
                             at their midpoint.
        n_points:            Number of sweep points per parameter (default: 50)
        threshold_variation: Fractional variation considered "saturated" (default: 0.05)
        threshold_fraction:  Fraction of range that must show saturation (default: 0.50)

    Returns:
        dict with keys:
            passed:     bool — True if no saturation detected in any parameter
            findings:   list of dicts, one per (parameter, output_key) pair:
                            parameter, output_key, saturated, fraction_saturated,
                            output_range_fraction, message
    """
    findings = []
    overall_passed = True

    # Build midpoint values for all parameters
    midpoints = {k: (v[0] + v[1]) / 2.0 for k, v in input_ranges.items()}

    for param, (lo, hi) in input_ranges.items():
        sweep_vals = np.linspace(lo, hi, n_points)
        outputs_by_key = {}

        for val in sweep_vals:
            kwargs = dict(midpoints)
            kwargs[param] = val
            result = func(**kwargs)

            if isinstance(result, dict):
                for key, out_val in result.items():
                    if np.isscalar(out_val) and np.isfinite(out_val):
                        outputs_by_key.setdefault(key, []).append(float(out_val))
            elif np.isscalar(result) and np.isfinite(result):
                outputs_by_key.setdefault("output", []).append(float(result))

        for key, vals in outputs_by_key.items():
            arr = np.array(vals)
            total_range = arr.max() - arr.min()

            if total_range == 0:
                frac_saturated = 1.0
            else:
                # Window spanning threshold_fraction of the sweep length.
                # A window is "flat" if its local variation is < threshold_variation
                # of the total output range. Count what fraction of windows are flat.
                window = max(2, int(threshold_fraction * n_points))
                n_windows = n_points - window + 1
                flat_count = sum(
                    1
                    for i in range(n_windows)
                    if (arr[i : i + window].max() - arr[i : i + window].min())
                    < threshold_variation * total_range
                )
                frac_saturated = flat_count / max(1, n_windows)

            # Saturated = there exists a contiguous window spanning >50% of the
            # input range over which the output varies by < 5% of total range.
            saturated = frac_saturated >= threshold_fraction
            if saturated:
                overall_passed = False

            msg = (
                f"SATURATED: '{key}' flat (< {threshold_variation*100:.0f}% total range) "
                f"in {frac_saturated*100:.0f}% of sliding windows "
                f"over {param} range [{lo:.3g}, {hi:.3g}]"
                if saturated
                else f"OK: '{key}' vs '{param}'"
            )

            findings.append(
                {
                    "parameter": param,
                    "output_key": key,
                    "saturated": saturated,
                    "fraction_saturated": frac_saturated,
                    "message": msg,
                }
            )

    return {"passed": overall_passed, "findings": findings}
