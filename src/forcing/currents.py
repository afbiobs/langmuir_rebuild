"""
Closed-basin surface current solver.

Solves the 1D vertically-resolved momentum balance for the wind-driven current
in a closed basin (no net throughflow). This is NOT the 3% rule; see AR-009.

The governing equation is:
    d/dz [nu_T(z) × du/dz] = 0    (z in [-h, 0])

with boundary conditions:
    nu_T × du/dz |_{z=0}   = tau_s / rho_water   (wind stress at surface)
    nu_T × du/dz |_{z=-h}  = 0                   (no bottom stress; stress-free)
    integral_{-h}^{0} u(z) dz = 0               (closed basin: zero net transport)

The stress-free bottom BC applies in the context of a closed-basin wind-driven
circulation where the bottom stress acts on the return flow at larger scale.
The zero-transport constraint imposes the pressure gradient that drives the
compensating return flow.

Source: AR-009 in assumptions_register.md.
"""

import numpy as np
from typing import Callable

# Physical constants
_RHO_WATER = 1000.0    # water density [kg/m³]
_RHO_AIR = 1.2         # air density [kg/m³]


def surface_velocity_1d(
    U10: float,
    depth: float,
    nu_T_profile: Callable[[np.ndarray], np.ndarray],
    drag_method: str = "coare35",
    n_grid: int = 200,
) -> float:
    """
    Near-surface velocity from 1D vertically-resolved momentum balance [m/s].

    Solves the wind-driven Couette-type problem in a closed basin using a
    cell-centred (staggered) grid to avoid the nu_T=0 boundary singularity.
    The returned value is the velocity at the uppermost cell centre (z = -dz/2),
    which approximates the surface velocity without the logarithmic singularity.

    Note: This function is provided for diagnostics. For the CL Rayleigh number,
    compute_forcing() uses the stress-based velocity scale (tau × h / nu_T_mean)
    which avoids the near-surface singularity entirely. See forcing/__init__.py.

    Parameters:
        U10:          Wind speed at 10 m height [m/s]
        depth:        Water depth h [m]
        nu_T_profile: Callable z -> nu_T [m²/s]
        drag_method:  "coare35" or "lake_low"
        n_grid:       Number of cell-centre grid points (default: 200)

    Returns:
        U_near_surface: Velocity at z = -dz/2 [m/s] (uppermost cell centre)

    Source: AR-009.
    """
    z, u = closed_basin_profile(U10, depth, nu_T_profile, drag_method, n_grid)
    return float(u[-1])  # uppermost cell centre


def closed_basin_profile(
    U10: float,
    depth: float,
    nu_T_profile: Callable[[np.ndarray], np.ndarray],
    drag_method: str = "coare35",
    n_grid: int = 200,
) -> tuple:
    """
    Full vertical current profile u(z) for the closed-basin 1D solution.

    Uses a cell-centred (staggered) grid: z_i = -h + (i+0.5) × dz for
    i = 0, ..., n_grid-1. This avoids evaluating nu_T at exactly z = 0 or
    z = -h where the parabolic profile is zero (logarithmic singularity in
    the velocity integral). nu_T at cell centres is always > 0.

    Parameters:
        U10:          Wind speed at 10 m height [m/s]
        depth:        Water depth h [m]
        nu_T_profile: Callable z -> nu_T [m²/s]
        drag_method:  "coare35" or "lake_low"
        n_grid:       Number of cell-centre points (default: 200)

    Returns:
        z:    Cell-centre coordinates [m], from -h+dz/2 to -dz/2
        u:    Current profile [m/s], satisfying integral(u dz) ≈ 0

    Notes:
        Verify zero-transport: trapz(u, z) should be near zero.
        The zero-transport constraint is approximate (cell-centred quadrature)
        but error is O(dz^2).

    Source: AR-009.
    """
    from .wind import drag_coefficient

    if U10 < 0:
        raise ValueError(f"U10 must be non-negative, got {U10:.3f} m/s")

    cd = drag_coefficient(U10, method=drag_method)
    tau_s = _RHO_AIR * cd * U10 ** 2
    bc_surface = tau_s / _RHO_WATER

    # Cell-centred grid: avoids nu_T = 0 at exact boundaries
    dz = depth / n_grid
    z = np.linspace(-depth + dz / 2, -dz / 2, n_grid)

    nu_T = nu_T_profile(z)
    nu_T = np.maximum(nu_T, 1e-9)  # molecular viscosity floor [m²/s]

    # Cumulative integral of 1/nu_T from bottom cell centre to each cell centre
    inv_nu = 1.0 / nu_T
    integral_inv_nu = np.zeros(n_grid)
    for i in range(1, n_grid):
        integral_inv_nu[i] = integral_inv_nu[i - 1] + 0.5 * (inv_nu[i - 1] + inv_nu[i]) * dz

    # Zero-transport constraint: C = -bc_surface × mean(integral_inv_nu)
    mean_integral = float(np.mean(integral_inv_nu))  # cell-centred average
    C = -bc_surface * mean_integral

    u = C + bc_surface * integral_inv_nu
    return z, u
