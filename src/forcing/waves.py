"""
Fetch-limited wave spectrum and Stokes drift profiles.

Implements JONSWAP fetch-limited relationships with shallow-water corrections,
plus two Stokes drift profile approximations (monochromatic and Webb & Fox-Kemper).

The differential Stokes drift D'(z) = du_s/dz [1/s] is the forcing input to the
CL instability solver.

IMPORTANT — Stokes drift formula:
    The correct formula for surface Stokes drift from a monochromatic wave is:
        u_s(0) = (H_s/2)^2 × omega_p × k_p   [m/s]
    This is derived from the first-order Stokes expansion of the wave velocity field.
    See failure_log.md FL-006 for the dimensional error in the implementation spec.

References:
    JONSWAP: Hasselmann et al. (1973), with fetch-limited forms from
             Young (1999) and Kahma & Calkoen (1992).
    Webb & Fox-Kemper: Webb, A. & Fox-Kemper, B. (2011). Ocean Modelling 37(1-2), 85-98.
    AR-012 in assumptions_register.md.
"""

import numpy as np
from typing import Optional

_G = 9.81        # gravitational acceleration [m/s²]
_GAMMA_J = 3.3   # JONSWAP peak enhancement factor [-]


def jonswap_parameters(U10: float, fetch: float, depth: float) -> dict:
    """
    Fetch-limited JONSWAP wave parameters with shallow-water correction.

    Uses the JONSWAP fetch-limited formulae (Young 1999):
        H_s = 1.6e-3 × (U10^2 / g) × X_tilde^0.5
        f_p = 3.5 × (g / U10) × X_tilde^{-1/3}
    where X_tilde = g × fetch / U10^2 is the dimensionless fetch.

    Peak wavenumber k_p is found by solving the dispersion relation:
        omega_p^2 = g × k_p × tanh(k_p × depth)
    using Newton iteration for finite-depth correction.

    Parameters:
        U10:   Wind speed at 10 m height [m/s]
        fetch: Wind fetch [m]
        depth: Water depth [m]

    Returns:
        dict with keys:
            H_s:       Significant wave height [m]
            T_p:       Peak wave period [s]
            f_p:       Peak wave frequency [Hz]
            omega_p:   Peak angular frequency [rad/s]
            k_p:       Peak wavenumber [1/m]
            lambda_p:  Peak wavelength [m]
            X_tilde:   Dimensionless fetch [-]
            shallow_water: bool, True if k_p × depth < pi (shallow wave condition)

    Source: AR-012, Young (1999).
    """
    if U10 <= 0:
        raise ValueError(f"U10 must be positive, got {U10:.3f} m/s")
    if fetch <= 0:
        raise ValueError(f"fetch must be positive, got {fetch:.1f} m")
    if depth <= 0:
        raise ValueError(f"depth must be positive, got {depth:.3f} m")

    X_tilde = _G * fetch / U10 ** 2  # dimensionless fetch [-]

    # Significant wave height [m]
    H_s = 1.6e-3 * (U10 ** 2 / _G) * X_tilde ** 0.5

    # Peak frequency [Hz]
    f_p = 3.5 * (_G / U10) * X_tilde ** (-1.0 / 3.0)
    omega_p = 2.0 * np.pi * f_p  # [rad/s]
    T_p = 1.0 / f_p               # [s]

    # Dispersion relation: omega^2 = g k tanh(k h)
    # Newton iteration starting from deep-water k = omega^2 / g
    k_p = omega_p ** 2 / _G  # deep-water initial guess
    for _ in range(30):
        kh = k_p * depth
        f = _G * k_p * np.tanh(kh) - omega_p ** 2
        df = _G * (np.tanh(kh) + kh / np.cosh(kh) ** 2)
        k_new = k_p - f / df
        if abs(k_new - k_p) < 1e-10:
            break
        k_p = k_new
    k_p = float(k_p)

    lambda_p = 2.0 * np.pi / k_p  # [m]
    shallow_water = bool(k_p * depth < np.pi)

    return {
        "H_s": float(H_s),
        "T_p": float(T_p),
        "f_p": float(f_p),
        "omega_p": float(omega_p),
        "k_p": float(k_p),
        "lambda_p": float(lambda_p),
        "X_tilde": float(X_tilde),
        "shallow_water": shallow_water,
    }


def stokes_drift_surface(U10: float, fetch: float, depth: float) -> float:
    """
    Surface Stokes drift from monochromatic (peak-frequency) approximation [m/s].

    u_s(0) = (H_s / 2)^2 × omega_p × k_p

    This is the standard linear-wave result for the surface Stokes drift of a
    monochromatic wave with amplitude a = H_s / 2. Units are [m/s]:
        [m]^2 × [rad/s] × [1/m] = [m/s].

    Parameters:
        U10:   Wind speed at 10 m height [m/s]
        fetch: Wind fetch [m]
        depth: Water depth [m]

    Returns:
        u_s0: Surface Stokes drift [m/s]

    Notes:
        This formula is NOT the same as (2π/T_p) × (π H_s / lambda_p), which has
        units [rad/s] not [m/s]. See failure_log.md FL-006.

    Source: Standard linear wave theory; Kenyon (1969).
    """
    p = jonswap_parameters(U10, fetch, depth)
    a = p["H_s"] / 2.0           # wave amplitude [m]
    u_s0 = a ** 2 * p["omega_p"] * p["k_p"]  # [m/s]
    return float(u_s0)


def stokes_drift_profile(
    z: np.ndarray,
    U10: float,
    fetch: float,
    depth: float,
    method: str = "webb_fox_kemper",
) -> np.ndarray:
    """
    Vertical profile of Stokes drift [m/s].

    Parameters:
        z:      Vertical coordinate [m], values in [-depth, 0]
        U10:    Wind speed at 10 m height [m/s]
        fetch:  Wind fetch [m]
        depth:  Water depth [m]
        method: "monochromatic" or "webb_fox_kemper" (default)

    Returns:
        u_s: Stokes drift profile [m/s], same shape as z

    Methods:
        "monochromatic":
            u_s(z) = u_s0 × exp(2 × k_p × z)
            WARNING: overestimates surface drift (ignores spectral broadening),
            underestimates drift at depth (ignores longer waves). Included as
            comparator only.

        "webb_fox_kemper":
            Exponential integral approximation for the JONSWAP spectrum:
                u_s(z) = u_s0 × exp(2 k_e z) / (1 − 8 k_e z)
            where k_e is an effective wavenumber defined such that the profile
            matches the JONSWAP spectral integral to leading order.
            For JONSWAP, k_e ≈ k_p (Webb & Fox-Kemper 2011, eq. 24).

    Source: Webb & Fox-Kemper (2011). Ocean Modelling 37(1-2), 85-98.
            AR-012.
    """
    z = np.asarray(z, dtype=float)
    p = jonswap_parameters(U10, fetch, depth)
    u_s0 = stokes_drift_surface(U10, fetch, depth)
    k_p = p["k_p"]

    if method == "monochromatic":
        return u_s0 * np.exp(2.0 * k_p * z)

    elif method == "webb_fox_kemper":
        # Effective wavenumber: k_e ≈ k_p for JONSWAP with gamma=3.3
        # The denominator (1 - 8 k_e z) is always > 1 for z < 0 (enhances depth penetration).
        k_e = k_p
        denom = 1.0 - 8.0 * k_e * z   # > 1 for z < 0
        # Guard against denominator becoming non-positive (deep z with low k_e)
        denom = np.maximum(denom, 0.1)
        return u_s0 * np.exp(2.0 * k_e * z) / denom

    else:
        raise ValueError(
            f"Unknown method: '{method}'. Use 'monochromatic' or 'webb_fox_kemper'."
        )


def differential_drift(
    z: np.ndarray,
    U10: float,
    fetch: float,
    depth: float,
    method: str = "webb_fox_kemper",
) -> np.ndarray:
    """
    Vertical derivative of Stokes drift: D'(z) = du_s/dz [1/s].

    This is the forcing input to the CL instability equations. Positive D'(z)
    means Stokes drift decreasing upward (as expected: drift is maximum at the
    surface).

    Computed by central differences on a fine grid, then interpolated to z.
    This avoids analytically differentiating the Webb & Fox-Kemper denominator.

    Parameters:
        z:      Vertical coordinate [m], values in [-depth, 0]
        U10:    Wind speed at 10 m height [m/s]
        fetch:  Wind fetch [m]
        depth:  Water depth [m]
        method: "monochromatic" or "webb_fox_kemper" (default)

    Returns:
        D_prime: du_s/dz [1/s], same shape as z. Positive values are the
                 physically expected sign: drift is maximum at the surface and
                 decays toward the bottom, so du_s/dz > 0 in the convention
                 where z=0 at surface and z=-h at bottom.

    Source: Derived from stokes_drift_profile; no independent reference needed.
    """
    z = np.asarray(z, dtype=float)

    # Fine grid for differentiation
    z_fine = np.linspace(-depth, 0.0, 500)
    dz = z_fine[1] - z_fine[0]
    u_s_fine = stokes_drift_profile(z_fine, U10, fetch, depth, method=method)

    # Central differences (interior), one-sided at endpoints
    D_prime_fine = np.gradient(u_s_fine, z_fine)

    # Interpolate to requested z points
    D_prime = np.interp(z, z_fine, D_prime_fine)
    return D_prime
