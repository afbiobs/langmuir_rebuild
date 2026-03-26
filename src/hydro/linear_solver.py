"""
Linear CL perturbation solver: small-l asymptotic expansion.

Implements H&P (2017) §3 to compute R₀, R*₂, and the linear neutral curve
R_L(l) = R₀ + l²R*₂ + (γR₀)/l² to leading order in l and γ.

Key results used:
  - R₀⁻¹ = ∫₋₁⁰ ψ̃₁ U' dz                        (eq. 26)
  - R*₂   = R₀ × ∫₋₁⁰ ψ̂₃ U' dz                   (γ=0 part of R₂)
  - R̃₂   = R₀  (linear case; ψ̃₃ = ψ̃₁ always)    (eq. 67)
  - l_cL  = (γ R₀ / |R*₂|)^(1/4)                   (eq. 67)
  - R_cL  = R₀ + 2(γ R₀ |R*₂|)^(1/2)               (eq. 70)

ψ̃₃ = ψ̃₁ always because both satisfy ψ̃'''' = D' with identical BCs.
This collapses R̃₂ = R₀ exactly (linear case), confirming eq. (67).

Reference: Hayes & Phillips (2017) §3; Cox & Leibovich (1993).
"""

import numpy as np
from dataclasses import dataclass
from scipy.optimize import minimize_scalar

from src.hydro.profiles import PolynomialProfile
from src.hydro.robin_bc import RobinBC


@dataclass(frozen=True)
class CriticalResult:
    """
    Critical point of a CL neutral curve.

    Fields:
        R_c:    Critical Rayleigh number [-]
        l_c:    Critical spanwise wavenumber [-]
        method: Solver provenance string
    """
    R_c: float   # [-]
    l_c: float   # [-]
    method: str


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _biharmonic_neumann(rhs: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    Solve ψ'''' = rhs on z ∈ [-1, 0] with ψ'' = ψ = 0 at both boundaries.

    Integrates four times using cumulative trapezoid, then applies 4 BCs
    to determine the 4 constants A, B, C, D in ψ = PI + A + Bz + Cz² + Dz³.

    Parameters:
        rhs: Right-hand side values on z grid [-]
        z:   Uniform grid on [-1, 0] [-]

    Returns:
        psi: Solution values on z grid [-]
    """
    n = len(z)
    dz = z[1] - z[0]

    # Build four cumulative integrals of rhs from z[0] = -1
    I = rhs.copy()
    for _ in range(4):
        out = np.zeros(n)
        for i in range(1, n):
            out[i] = out[i - 1] + 0.5 * (I[i] + I[i - 1]) * dz
        I = out
    I4 = I  # particular integral

    # Recover I2 (second integral of rhs) for ψ'' BC
    I2 = rhs.copy()
    for _ in range(2):
        out = np.zeros(n)
        for i in range(1, n):
            out[i] = out[i - 1] + 0.5 * (I2[i] + I2[i - 1]) * dz
        I2 = out

    # BCs: ψ(0)=0, ψ(-1)=0, ψ''(0)=0, ψ''(-1)=0
    # ψ  = I4 + A + Bz + Cz² + Dz³
    # ψ'' = I2 + 2C + 6Dz
    #
    # ψ(0)=0   → I4[-1] + A = 0               → A = -I4[-1]
    # ψ''(0)=0 → I2[-1] + 2C = 0              → C = -I2[-1]/2
    # ψ''(-1)=0→ I2[0] + 2C - 6D = 0 (I2[0]=0)→ D = 2C/6 = C/3
    # ψ(-1)=0  → I4[0] + A - B + C - D = 0 (I4[0]=0)
    #           → A - B + C - D = 0            → B = A + C - D

    A = -I4[-1]
    C = -I2[-1] / 2.0
    D = C / 3.0
    B = A + C - D

    return I4 + A + B * z + C * z ** 2 + D * z ** 3


def _d2(f: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Second derivative by central differences; one-sided at boundaries."""
    dz = z[1] - z[0]
    d2 = np.empty_like(f)
    d2[1:-1] = (f[2:] - 2 * f[1:-1] + f[:-2]) / dz ** 2
    d2[0] = d2[1]
    d2[-1] = d2[-2]
    return d2


def _cumtrapz(f: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Cumulative trapezoid integral from z[0]."""
    n = len(z)
    out = np.zeros(n)
    for i in range(1, n):
        out[i] = out[i - 1] + 0.5 * (f[i] + f[i - 1]) * (z[i] - z[i - 1])
    return out


def _precompute(
    D_prime: PolynomialProfile,
    U_prime: PolynomialProfile,
    n_pts: int = 2000,
) -> dict:
    """
    Compute R₀ and R*₂ for given profiles (γ-independent quantities).

    Returns dict with keys: z, R0, R_star_2, psi_tilde_1
    """
    z = np.linspace(-1.0, 0.0, n_pts)
    D_vals = D_prime(z)
    U_vals = U_prime(z)

    # Step 1: ψ̃₁'''' = D', Neumann BCs
    psi_tilde_1 = _biharmonic_neumann(D_vals, z)

    # Step 2: R₀ from eq. (26)
    R0 = 1.0 / np.trapz(psi_tilde_1 * U_vals, z)

    # Step 3: ψ₁ = -R₀ ψ̃₁
    psi_1 = -R0 * psi_tilde_1

    # Step 4: u₂ (γ=0; c₀=0 by solvability, c₁ = -∫cum2 dz)
    # u₂'' = u₀ + U'ψ₁  → integrate twice, subtract mean
    integrand = np.ones(n_pts) + U_vals * psi_1
    cum1 = _cumtrapz(integrand, z)
    cum2 = _cumtrapz(cum1, z)
    u_2 = cum2 - np.trapz(cum2, z)  # enforce ∫u₂ dz = 0

    # Step 5: ψ̂₃'''' = 2ψ₁'' - R₀D'u₂, Neumann BCs
    rhs_hat3 = 2.0 * _d2(psi_1, z) - R0 * D_vals * u_2
    psi_hat_3 = _biharmonic_neumann(rhs_hat3, z)

    # Step 6: R*₂ = R₀ × ∫ψ̂₃ U' dz
    R_star_2 = R0 * np.trapz(psi_hat_3 * U_vals, z)

    return {"z": z, "R0": R0, "R_star_2": R_star_2, "psi_tilde_1": psi_tilde_1}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_R0(
    D_prime: PolynomialProfile,
    U_prime: PolynomialProfile,
    n_pts: int = 2000,
) -> float:
    """
    Leading Rayleigh coefficient R₀ [-].

    R₀⁻¹ = ∫₋₁⁰ ψ̃₁(z) U'(z) dz   (H&P 2017 eq. 26)

    For uniform profiles D' = U' = 1: R₀ = 120 exactly.

    Parameters:
        D_prime: Differential drift profile D'(z) [-]
        U_prime: Mean shear profile U'(z) [-]
        n_pts:   Grid resolution for numerical integration

    Returns:
        R0: Leading Rayleigh coefficient [-]

    Source: H&P (2017) eq. (26).
    """
    return _precompute(D_prime, U_prime, n_pts)["R0"]


def neutral_curve_linear(
    l_array: np.ndarray,
    D_prime: PolynomialProfile,
    U_prime: PolynomialProfile,
    bc: RobinBC,
    n_pts: int = 2000,
) -> np.ndarray:
    """
    Linear neutral curve R_L(l) from the leading-order small-l expansion.

    R_L(l) = R₀ + l² R*₂ + (γ R₀) / l²

    Valid near l_cL = O(γ^{1/4}). For Neumann BCs (γ=0), R_L(l) = R₀ + l²R*₂.

    Parameters:
        l_array: Spanwise wavenumber values [-]
        D_prime: Differential drift profile D'(z) [-]
        U_prime: Mean shear profile U'(z) [-]
        bc:      Robin boundary conditions [-]
        n_pts:   Grid resolution

    Returns:
        R_linear: Neutral curve values [-]

    Source: H&P (2017) eqs. (5b), (65); γ̃ scaling from eq. (6).
    """
    l_array = np.asarray(l_array, dtype=float)
    p = _precompute(D_prime, U_prime, n_pts)
    R0, R_star_2 = p["R0"], p["R_star_2"]
    gamma = bc.gamma

    R = R0 + l_array ** 2 * R_star_2
    if gamma > 0.0:
        R = R + gamma * R0 / l_array ** 2
    return R


def critical_linear(
    D_prime: PolynomialProfile,
    U_prime: PolynomialProfile,
    bc: RobinBC,
    n_pts: int = 2000,
    l_search_range: tuple = (1e-4, 3.0),
) -> CriticalResult:
    """
    Critical linear Rayleigh number R_cL and wavenumber l_cL.

    For Robin BCs (γ > 0), uses the closed-form minimum of R_L(l):
        l_cL  = (γ R₀ / |R*₂|)^{1/4}        (H&P 2017 eq. 67)
        R_cL  = R₀ + 2(γ R₀ |R*₂|)^{1/2}    (H&P 2017 eq. 70)

    For Neumann BCs (γ=0), l_cL = 0 and R_cL = R₀.

    Parameters:
        D_prime, U_prime: Profile functions [-]
        bc:               Robin boundary conditions [-]
        n_pts:            Grid resolution
        l_search_range:   Fallback numerical search range (not used for γ>0)

    Returns:
        CriticalResult with R_c = R_cL, l_c = l_cL

    Source: H&P (2017) eqs. (67), (70).
    """
    p = _precompute(D_prime, U_prime, n_pts)
    R0, R_star_2 = p["R0"], p["R_star_2"]
    gamma = bc.gamma

    if bc.is_neumann or gamma == 0.0:
        return CriticalResult(R_c=R0, l_c=0.0, method="asymptotic_linear_neumann")

    # Closed-form critical point of R₀ + l²R*₂ + γR₀/l²
    # dR/dl = 2l R*₂ - 2γR₀/l³ = 0  →  l⁴ = γR₀/|R*₂|
    abs_R_star_2 = abs(R_star_2)
    l_c = (gamma * R0 / abs_R_star_2) ** 0.25
    R_c = R0 + 2.0 * (gamma * R0 * abs_R_star_2) ** 0.5

    return CriticalResult(R_c=R_c, l_c=l_c, method="asymptotic_linear")
