"""
Forcing module: wind, waves, currents, and eddy viscosity.

Entry point: compute_forcing(U10, depth, fetch, ...) -> ForcingState

All quantities in SI units. The ForcingState frozen dataclass captures the
complete forcing description at a single instant, with full provenance metadata.

Key physics decisions encoded here:
- Water-side u*_water is used in nu_T and Ra (not air-side u*_air)
- Closed-basin surface current (not the 3% oceanic rule)
- Correct Stokes drift formula: u_s(0) = (H_s/2)^2 * omega_p * k_p [m/s]
- La_t = sqrt(u*_water / U_s0) (McWilliams 1997 water-side convention)

See assumptions_register.md AR-009 through AR-012 for sourcing of each choice.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional

from .wind import friction_velocity, drag_coefficient
from .waves import jonswap_parameters, stokes_drift_surface, stokes_drift_profile, differential_drift
from .eddy_viscosity import parabolic_nu_T, representative_nu_T
from .currents import surface_velocity_1d

# Density ratio for air-to-water friction velocity conversion
_RHO_AIR = 1.2          # [kg/m³]
_RHO_WATER = 1000.0     # [kg/m³]
_DENSITY_RATIO = np.sqrt(_RHO_AIR / _RHO_WATER)   # u*_water / u*_air ≈ 0.0346


@dataclass(frozen=True)
class ForcingState:
    """
    Complete forcing description at a single time. All dimensional, SI units.

    This is the output of compute_forcing() and the input to the hydro module.
    Every field has units documented. No field may be NaN after a successful
    compute_forcing() call.

    Provenance fields (drag_method, drift_method) record which parameterisation
    was used, enabling reproducibility audits.
    """
    U10: float                              # Wind speed at 10 m height [m/s]
    u_star_air: float                       # Atmospheric friction velocity [m/s]
    u_star_water: float                     # Water-side friction velocity [m/s]
    U_surface: float                        # Surface current from closed-basin solver [m/s]
    stokes_drift_surface: float             # Surface Stokes drift u_s(0) [m/s]
    stokes_drift_profile: np.ndarray        # Stokes drift profile u_s(z) [m/s] on z_grid
    differential_drift_profile: np.ndarray  # D'(z) = du_s/dz [1/s] on z_grid
    z_grid: np.ndarray                      # Vertical coordinate [m], from -depth to 0
    depth: float                            # Water depth h [m]
    fetch: float                            # Wind fetch [m]
    H_s: float                              # Significant wave height [m]
    T_p: float                              # Peak wave period [s]
    f_p: float                              # Peak wave frequency [Hz]
    omega_p: float                          # Peak angular frequency [rad/s]
    k_p: float                              # Peak wavenumber [1/m]
    lambda_p: float                         # Peak wavelength [m]
    X_tilde: float                          # Dimensionless fetch [-]
    wave_steepness: float                   # ak = (H_s / 2) × k_p [-]
    La_t: float                             # Turbulent Langmuir number (water-side) [-]
    nu_T: float                             # Depth-averaged eddy viscosity [m²/s]
    nu_T_profile: np.ndarray               # nu_T(z) profile [m²/s] on z_grid
    Ra: float                               # Langmuir–Rayleigh number [-]
    timestamp: datetime                     # Observation time (UTC)
    drag_method: str                        # Provenance: drag parameterisation used
    drift_method: str                       # Provenance: Stokes drift method used

    def __post_init__(self):
        # Validate no NaN in scalar fields
        scalar_fields = [
            "U10", "u_star_air", "u_star_water", "U_surface",
            "stokes_drift_surface", "depth", "fetch", "H_s", "T_p",
            "f_p", "omega_p", "k_p", "lambda_p", "X_tilde", "wave_steepness",
            "La_t", "nu_T", "Ra",
        ]
        for field in scalar_fields:
            val = getattr(self, field)
            if not np.isfinite(val):
                raise ValueError(
                    f"ForcingState.{field} = {val} (must be finite). "
                    "Check forcing inputs."
                )


def compute_forcing(
    U10: float,
    depth: float,
    fetch: float,
    timestamp: Optional[datetime] = None,
    drag_method: str = "coare35",
    drift_method: str = "webb_fox_kemper",
    n_grid: int = 200,
) -> ForcingState:
    """
    Compute the complete forcing state for given wind and site conditions.

    This is the primary entry point to the forcing module. It orchestrates
    wind, waves, currents, and eddy viscosity computations and assembles them
    into a ForcingState.

    Parameters:
        U10:         Wind speed at 10 m height [m/s]
        depth:       Water depth h [m]
        fetch:       Wind fetch F [m]
        timestamp:   Observation datetime (UTC); defaults to now if None
        drag_method: "coare35" (default) or "lake_low"
        drift_method: "webb_fox_kemper" (default) or "monochromatic"
        n_grid:      Vertical grid resolution (default: 200)

    Returns:
        ForcingState with all fields populated.

    Physics chain:
        1. u*_air = sqrt(C_D) × U10                     [wind.py]
        2. u*_water = sqrt(rho_air/rho_water) × u*_air  [density ratio]
        3. nu_T(z) = kappa × u*_water × (z+h)(1-(z+h)/h) [eddy_viscosity.py]
        4. nu_T_mean = integral(nu_T) / h                [eddy_viscosity.py]
        5. U_surface from closed-basin 1D solver         [currents.py]
        6. H_s, T_p, k_p from JONSWAP(U10, fetch, depth) [waves.py]
        7. H_s, T_p, f_p, omega_p, k_p, lambda_p         [waves.py]
        8. u_s(0) = (H_s/2)^2 × omega_p × k_p           [waves.py]
        9. u_s(z) from drift_method                      [waves.py]
       10. D'(z) = d u_s/dz                              [waves.py]
       11. La_t = sqrt(u*_water / u_s(0))                [McWilliams 1997]
       12. Ra = U_surface × u_s(0) × h^2 / nu_T_mean^2   [rayleigh.py]

    Notes:
        Ra here uses U_surface as the wind-driven shear velocity scale and
        u_s(0) as the Stokes drift scale. The Rayleigh number is dimensionless
        by construction. Typical values for Lough Neagh (h=9m) are O(10,000)
        at moderate wind, deeply supercritical (R_cNL ≈ 122).

    Source: AR-009 through AR-012.
    """
    if timestamp is None:
        timestamp = datetime.now(timezone.utc)

    # Step 1-2: friction velocities
    u_star_air = friction_velocity(U10, method=drag_method)
    u_star_water = _DENSITY_RATIO * u_star_air

    # Step 3-4: eddy viscosity
    z_grid = np.linspace(-depth, 0.0, n_grid)
    nu_T_prof = parabolic_nu_T(z_grid, u_star_water, depth)
    nu_T_mean = representative_nu_T(nu_T_prof, z_grid)

    # Guard against degenerate nu_T (e.g. U10 = 0)
    if nu_T_mean <= 0.0:
        nu_T_mean = 1e-6

    # Step 5: wind-driven velocity scale (stress-based).
    # The parabolic nu_T profile has nu_T -> 0 at both boundaries, which makes
    # the 1D integral int(dz/nu_T) diverge and the surface velocity ill-defined
    # (see failure_log.md FL-007 for the singularity analysis). The physical
    # velocity scale entering the CL Rayleigh number is the SHEAR amplitude:
    #
    #   U_surface = tau_wind × h / (rho_water × nu_T_mean)
    #             = bc_surface × h / nu_T_mean   [m/s]
    #
    # This is the velocity scale of the depth-varying wind-driven current in the
    # bulk of the water column, avoiding the surface singularity. It equals the
    # exact surface velocity for a Couette flow with constant viscosity nu_T_mean,
    # and is the natural scaling for the CL shear term.
    #
    # Source: Momentum balance in closed basin; see AR-009.
    _RHO_AIR_LOC = 1.2
    _RHO_WATER_LOC = 1000.0
    cd_loc = drag_coefficient(U10, method=drag_method)
    bc_surface = _RHO_AIR_LOC * cd_loc * U10 ** 2 / _RHO_WATER_LOC  # [m^2/s^2]
    U_surface = bc_surface * depth / nu_T_mean  # [m/s]

    # Step 6-9: wave parameters and Stokes drift
    wave_params = jonswap_parameters(U10, fetch, depth)
    H_s = wave_params["H_s"]
    T_p = wave_params["T_p"]
    f_p = wave_params["f_p"]
    omega_p = wave_params["omega_p"]
    k_p = wave_params["k_p"]
    lambda_p = wave_params["lambda_p"]
    X_tilde = wave_params["X_tilde"]
    wave_steepness = 0.5 * H_s * k_p

    u_s0 = stokes_drift_surface(U10, fetch, depth)
    u_s_prof = stokes_drift_profile(z_grid, U10, fetch, depth, method=drift_method)
    D_prime = differential_drift(z_grid, U10, fetch, depth, method=drift_method)

    # Step 10: turbulent Langmuir number (water-side)
    # La_t = sqrt(u*_water / u_s0); guard against u_s0 = 0
    if u_s0 > 0:
        La_t = np.sqrt(u_star_water / u_s0)
    else:
        La_t = float("inf")

    # Step 11: Rayleigh number
    # Ra = U_surface × D_max × h^2 / nu_T_mean^2
    # D_max = max(u_s(z)) = u_s0 (surface value is maximum)
    D_max = u_s0
    Ra = abs(U_surface) * D_max * depth ** 2 / nu_T_mean ** 2

    return ForcingState(
        U10=float(U10),
        u_star_air=float(u_star_air),
        u_star_water=float(u_star_water),
        U_surface=float(U_surface),
        stokes_drift_surface=float(u_s0),
        stokes_drift_profile=u_s_prof,
        differential_drift_profile=D_prime,
        z_grid=z_grid,
        depth=float(depth),
        fetch=float(fetch),
        H_s=float(H_s),
        T_p=float(T_p),
        f_p=float(f_p),
        omega_p=float(omega_p),
        k_p=float(k_p),
        lambda_p=float(lambda_p),
        X_tilde=float(X_tilde),
        wave_steepness=float(wave_steepness),
        La_t=float(La_t),
        nu_T=float(nu_T_mean),
        nu_T_profile=nu_T_prof,
        Ra=float(Ra),
        timestamp=timestamp,
        drag_method=drag_method,
        drift_method=drift_method,
    )
