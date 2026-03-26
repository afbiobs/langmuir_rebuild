"""
Eddy viscosity profiles for the wind-driven surface layer.

Uses the parabolic mixing-length model (AR-010). All functions operate with
water-side friction velocity u*_water (not the air-side value).

Reference: Prandtl mixing-length theory; see assumptions_register.md AR-010.
"""

import numpy as np


def parabolic_nu_T(
    z: np.ndarray,
    u_star_water: float,
    depth: float,
    kappa: float = 0.41,
) -> np.ndarray:
    """
    Parabolic eddy viscosity profile [m²/s].

    ν_T(z) = κ × u*_water × (z + h) × (1 − (z + h) / h)

    where z ∈ [-h, 0], h = depth.

    This profile is zero at both boundaries (z = −h and z = 0) and peaks at
    mid-depth with maximum ν_T_max = κ × u*_water × h / 4.

    Parameters:
        z:            Vertical coordinate array [m], values in [-depth, 0]
        u_star_water: Water-side friction velocity [m/s]
        depth:        Water depth h [m], positive value
        kappa:        von Kármán constant [-] (default: 0.41)

    Returns:
        nu_T: Eddy viscosity profile [m²/s], same shape as z

    Notes:
        MUST use water-side u*_water, not the atmospheric u*_air.
        u*_water = sqrt(rho_air / rho_water) × u*_air ≈ 0.0346 × u*_air.
        Depth-average of ν_T(z) over [-h, 0] = κ × u*_water × h / 6.

    Source: AR-010, Prandtl mixing length.
    """
    z = np.asarray(z, dtype=float)
    if u_star_water < 0:
        raise ValueError(f"u_star_water must be non-negative, got {u_star_water:.6f}")
    if depth <= 0:
        raise ValueError(f"depth must be positive, got {depth:.3f}")

    zeta = z + depth          # zeta in [0, h]; 0 at bottom, h at surface
    nu_T = kappa * u_star_water * zeta * (1.0 - zeta / depth)

    # Enforce non-negative (floating-point clips at boundaries)
    nu_T = np.maximum(nu_T, 0.0)
    return nu_T


def representative_nu_T(profile: np.ndarray, z_grid: np.ndarray) -> float:
    """
    Depth-averaged eddy viscosity [m²/s] for use in Rayleigh number computation.

    Integrates nu_T(z) over the full depth using the trapezoidal rule.

    Parameters:
        profile: nu_T(z) values [m²/s] on z_grid
        z_grid:  Vertical coordinate array [m], values in [-depth, 0]

    Returns:
        nu_T_mean: Depth-averaged eddy viscosity [m²/s]

    Notes:
        For the parabolic profile, the exact depth-average is kappa * u_star_water * h / 6.
        This function integrates numerically to remain general.

    Source: AR-010.
    """
    profile = np.asarray(profile, dtype=float)
    z_grid = np.asarray(z_grid, dtype=float)
    depth = float(z_grid.max() - z_grid.min())
    if depth <= 0:
        raise ValueError("z_grid must span a positive depth range.")
    return float(np.trapz(profile, z_grid) / depth)
