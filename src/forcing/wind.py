"""
Wind stress and friction velocity.

Implements COARE 3.5 (Fairall et al. 2003) and a sheltered-lake alternative.

All functions use SI units throughout. Air-side friction velocity is computed
from the neutral drag coefficient. Water-side u* is derived in __init__.py via
the density ratio sqrt(rho_air / rho_water).

Source: Fairall, C.W. et al. (2003). J. Climate 16(4), 571–591.
        COARE 3.5 update: Edson et al. (2013). J. Phys. Oceanogr. 43(8), 1589–1610.
"""

import numpy as np

# Physical constants
_KAPPA_VON_KARMAN = 0.40      # von Kármán constant [-] (COARE convention: 0.40)
_G = 9.81                      # gravitational acceleration [m/s²]
_Z_REF = 10.0                  # reference height [m]
_NU_AIR = 1.5e-5               # kinematic viscosity of air [m²/s]

# COARE 3.5 Charnock parameters (Edson et al. 2013)
_CHARNOCK_SLOPE = 0.017        # m = 0.017 m⁻¹s (slope in α = m × U10N + b)
_CHARNOCK_INTERCEPT = -0.005   # b = -0.005 (intercept)
_CHARNOCK_MIN = 1.8e-4         # minimum roughness (smooth flow limit) [-]

# Lake-low drag coefficient (sheltered, U10 < 5 m/s)
_CD_LAKE_LOW = 1.0e-3          # [-]


def drag_coefficient(U10: float, method: str = "coare35") -> float:
    """
    Neutral 10-m drag coefficient C_DN [-].

    Parameters:
        U10:    Wind speed at 10 m height [m/s]
        method: "coare35" (default) or "lake_low"

    Returns:
        C_DN: Neutral drag coefficient [-]

    Methods:
        "coare35": COARE 3.5 iterative Charnock formulation.
            Roughness: z0 = (α_c / g) × u*² + 0.11 × ν_air / u*
            Charnock: α_c = max(m × U10N + b, α_min) where m=0.017, b=-0.005
            Neutral: C_DN = [κ / ln(z_ref / z0)]²
            Iterated to convergence (typically 5 iterations).
        "lake_low": Sheltered-lake regime.
            C_DN = 1.0e-3 for U10 < 5 m/s.
            Blends to COARE 3.5 for U10 ≥ 5 m/s.

    Source: Edson et al. (2013) for COARE 3.5 coefficients.
            AR-009 in assumptions_register.md.
    """
    if U10 < 0:
        raise ValueError(f"U10 must be non-negative, got {U10:.3f} m/s")

    if method == "lake_low":
        if U10 < 5.0:
            return _CD_LAKE_LOW
        # Blend: at U10=5 lake_low gives 1e-3; COARE can give higher values.
        # Use COARE above 5 m/s.
        return drag_coefficient(U10, method="coare35")

    if method != "coare35":
        raise ValueError(f"Unknown drag method: '{method}'. Use 'coare35' or 'lake_low'.")

    if U10 < 0.1:
        # Very low wind: use smooth-flow limit
        return _CD_LAKE_LOW

    # Iterative COARE 3.5
    u_star = 0.035 * U10  # initial guess
    for _ in range(10):
        alpha_c = max(_CHARNOCK_SLOPE * U10 + _CHARNOCK_INTERCEPT, _CHARNOCK_MIN)
        z0 = alpha_c * u_star ** 2 / _G + 0.11 * _NU_AIR / (u_star + 1e-10)
        z0 = max(z0, 1e-6)  # numerical floor
        cd = (_KAPPA_VON_KARMAN / np.log(_Z_REF / z0)) ** 2
        u_star_new = np.sqrt(cd) * U10
        if abs(u_star_new - u_star) < 1e-8:
            break
        u_star = u_star_new

    return float(cd)


def friction_velocity(U10: float, method: str = "coare35") -> float:
    """
    Atmospheric friction velocity u*_air [m/s].

    u*_air = sqrt(C_DN) × U10

    Parameters:
        U10:    Wind speed at 10 m height [m/s]
        method: "coare35" (default) or "lake_low"

    Returns:
        u_star_air: Atmospheric friction velocity [m/s]

    Notes:
        This is the AIR-SIDE friction velocity. The water-side u*_water is
        computed in forcing/__init__.py as:
            u*_water = sqrt(rho_air / rho_water) × u*_air
                     ≈ 0.0346 × u*_air
        for rho_air = 1.2 kg/m³, rho_water = 1000 kg/m³.

        The water-side u*_water is used in the eddy viscosity profile (AR-010)
        and in the Rayleigh number (AR-009). Using the air-side value in these
        expressions is a dimensional error.

    Source: Fairall et al. (2003), Edson et al. (2013).
    """
    if U10 < 0:
        raise ValueError(f"U10 must be non-negative, got {U10:.3f} m/s")
    cd = drag_coefficient(U10, method=method)
    return float(np.sqrt(cd) * U10)
