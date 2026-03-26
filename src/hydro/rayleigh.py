"""
Rayleigh number computation and regime classification for Langmuir circulation.

The Langmuir-Rayleigh number Ra controls the onset and intensity of LC:
    Ra = U_surface × D_max × h^2 / nu_T^2

where:
    U_surface: wind-driven surface current [m/s]
    D_max:     maximum Stokes drift (surface value) [m/s]
    h:         water depth [m]
    nu_T:      representative (depth-averaged) eddy viscosity [m^2/s]

Onset occurs at Ra = R_cNL ≈ 122 for uniform profiles with small Robin BCs.
For typical Lough Neagh conditions (h=9m, U10=5 m/s), Ra ~ 10,000–50,000,
deeply supercritical. The regime classification reflects this.

Reference: Hayes & Phillips (2017) §2–3; Cox & Leibovich (1993).
"""

import numpy as np
from typing import Callable, Optional


def compute_rayleigh(
    U_surface: float,
    D_max: float,
    depth: float,
    nu_T: float,
) -> float:
    """
    Langmuir–Rayleigh number Ra [-].

    Ra = U_surface × D_max × h^2 / nu_T^2

    Parameters:
        U_surface: Wind-driven surface current [m/s]
        D_max:     Maximum Stokes drift (surface value, u_s(0)) [m/s]
        depth:     Water depth h [m]
        nu_T:      Representative (depth-averaged) eddy viscosity [m^2/s]

    Returns:
        Ra: Rayleigh number [-]

    Notes:
        Ra is always non-negative (absolute value of U_surface is used).
        If nu_T = 0 (e.g. U10 = 0), Ra = 0.

    Source: Hayes & Phillips (2017) eq. (1), §2.1.
    """
    if nu_T <= 0:
        return 0.0
    return float(abs(U_surface) * abs(D_max) * depth ** 2 / nu_T ** 2)


def classify_regime(Ra: float, R0: float = 120.0, RcNL: float = 122.194) -> str:
    """
    Classify the LC dynamical regime from the Rayleigh number.

    Parameters:
        Ra:   Rayleigh number [-]
        R0:   Onset threshold (asymptotic, gamma -> 0) [-] (default: 120)
        RcNL: Nonlinear critical Rayleigh number [-] (default: 122.194)

    Returns:
        regime: One of "subcritical", "near_onset", "moderate", "supercritical"

    Regime definitions:
        "subcritical"   : Ra < R0         — no LC predicted
        "near_onset"    : R0 <= Ra < 1.5 × RcNL  — onset physics dominates;
                          spacing sensitive to Ra, near l_cNL
        "moderate"      : 1.5 × RcNL <= Ra < 5 × RcNL  — finite-amplitude cells;
                          spacing weakly depends on Ra; coarsening likely active
        "supercritical" : Ra >= 5 × RcNL  — strongly forced; aspect ratio
                          approaching cap; coarsening maximised

    Notes:
        For Lough Neagh at U10 = 5 m/s, Ra ~ 20,000 >> 5 × RcNL ≈ 611.
        Nearly all realistic observations are in the "supercritical" regime.
        This is physically correct and does NOT indicate model saturation.

    Source: Regime labels are new to this model; thresholds from H&P (2017) §8.
    """
    if Ra < R0:
        return "subcritical"
    elif Ra < 1.5 * RcNL:
        return "near_onset"
    elif Ra < 5.0 * RcNL:
        return "moderate"
    else:
        return "supercritical"


def unstable_band(
    Ra: float,
    neutral_curve: Callable[[np.ndarray], np.ndarray],
    l_array: np.ndarray,
) -> tuple:
    """
    Unstable wavenumber band [l_min, l_max] where Ra > R_bar(l).

    For Ra below the neutral curve everywhere, returns (nan, nan).

    Parameters:
        Ra:            Rayleigh number [-]
        neutral_curve: Callable l_array -> R_bar(l), the nonlinear neutral curve [-]
        l_array:       Wavenumber array to evaluate on [normalised, -]

    Returns:
        (l_min, l_max): bounds of unstable band [normalised wavenumber, -]
                        Both nan if Ra is subcritical everywhere.

    Notes:
        The unstable band is where Ra exceeds the neutral curve R_bar(l).
        Within this band, finite-amplitude LC cells are sustained at those
        wavenumbers (hence widths L = 2pi/l).

    Source: H&P (2017) Figure 4.
    """
    l_array = np.asarray(l_array, dtype=float)
    R_neutral = neutral_curve(l_array)
    unstable = Ra > R_neutral

    if not np.any(unstable):
        return (float("nan"), float("nan"))

    l_unstable = l_array[unstable]
    return (float(l_unstable.min()), float(l_unstable.max()))
