"""
Simple baseline predictors for windrow spacing.

These baselines define the performance floor that physics-based candidates must beat.
They require no physics — only the observation statistics and OLS regression.

All functions return a dict containing:
    predicted:   np.ndarray of predicted spacings [m]
    coefficients: dict of fitted parameters with units
    rmse:        float [m]
    mae:         float [m]
"""

import numpy as np
from typing import Optional


def baseline_constant(observations: np.ndarray) -> dict:
    """
    Predict the mean observed spacing for all cases.

    This is the null model: it encodes only the marginal distribution of the
    target variable. Any physics-based model must beat this.

    Parameters:
        observations: Observed windrow spacings [m]

    Returns:
        dict with keys:
            predicted:    array of mean spacing repeated N times [m]
            coefficients: {"mean_spacing": float [m]}
            rmse:         float [m]
            mae:          float [m]
    """
    obs = np.asarray(observations, dtype=float)
    mean_spacing = float(np.mean(obs))
    predicted = np.full(len(obs), mean_spacing)

    residuals = predicted - obs
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    mae = float(np.mean(np.abs(residuals)))

    return {
        "predicted": predicted,
        "coefficients": {"mean_spacing_m": mean_spacing},
        "rmse": rmse,
        "mae": mae,
    }


def baseline_linear_wind(
    wind_speeds: np.ndarray,
    observations: np.ndarray,
) -> dict:
    """
    OLS regression of spacing on wind speed: spacing = a + b × U10.

    This captures any linear trend between wind and spacing. A physics-based
    model must produce a better fit than this two-parameter baseline.

    Parameters:
        wind_speeds:  Representative U10 for each observation [m/s]
        observations: Observed windrow spacings [m]

    Returns:
        dict with keys:
            predicted:    OLS predictions [m]
            coefficients: {"intercept_m": float, "slope_m_per_ms": float}
            rmse:         float [m]
            mae:          float [m]
            r_squared:    float [-]
    """
    U = np.asarray(wind_speeds, dtype=float)
    obs = np.asarray(observations, dtype=float)

    if len(U) != len(obs):
        raise ValueError(
            f"Length mismatch: wind_speeds={len(U)}, observations={len(obs)}"
        )

    # OLS: [intercept, slope] via normal equations
    A = np.column_stack([np.ones(len(U)), U])
    coeffs, _, _, _ = np.linalg.lstsq(A, obs, rcond=None)
    intercept, slope = float(coeffs[0]), float(coeffs[1])

    predicted = intercept + slope * U
    residuals = predicted - obs
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((obs - obs.mean()) ** 2))
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "predicted": predicted,
        "coefficients": {
            "intercept_m": intercept,
            "slope_m_per_ms": slope,
        },
        "rmse": float(np.sqrt(np.mean(residuals ** 2))),
        "mae": float(np.mean(np.abs(residuals))),
        "r_squared": r_sq,
    }


def baseline_depth_scaled(
    wind_speeds: np.ndarray,
    depths: np.ndarray,
    observations: np.ndarray,
) -> dict:
    """
    OLS regression of spacing on wind speed and depth:
        spacing = a + b × U10 + c × depth

    This is the minimal physics-aware baseline: it encodes the empirical fact
    that LC cell spacing scales with depth. A physics-based model must produce
    better predictions than this three-parameter baseline.

    Parameters:
        wind_speeds:  Representative U10 for each observation [m/s]
        depths:       Water depth at each observation site [m]
        observations: Observed windrow spacings [m]

    Returns:
        dict with keys:
            predicted:    OLS predictions [m]
            coefficients: {"intercept_m": float, "slope_wind_m_per_ms": float,
                           "slope_depth_m_per_m": float}
            rmse:         float [m]
            mae:          float [m]
            r_squared:    float [-]
    """
    U = np.asarray(wind_speeds, dtype=float)
    h = np.asarray(depths, dtype=float)
    obs = np.asarray(observations, dtype=float)

    if not (len(U) == len(h) == len(obs)):
        raise ValueError(
            f"Length mismatch: wind_speeds={len(U)}, depths={len(h)}, "
            f"observations={len(obs)}"
        )

    A = np.column_stack([np.ones(len(U)), U, h])
    coeffs, _, _, _ = np.linalg.lstsq(A, obs, rcond=None)
    intercept = float(coeffs[0])
    slope_wind = float(coeffs[1])
    slope_depth = float(coeffs[2])

    predicted = intercept + slope_wind * U + slope_depth * h
    residuals = predicted - obs
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((obs - obs.mean()) ** 2))
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "predicted": predicted,
        "coefficients": {
            "intercept_m": intercept,
            "slope_wind_m_per_ms": slope_wind,
            "slope_depth_m_per_m": slope_depth,
        },
        "rmse": float(np.sqrt(np.mean(residuals ** 2))),
        "mae": float(np.mean(np.abs(residuals))),
        "r_squared": r_sq,
    }
