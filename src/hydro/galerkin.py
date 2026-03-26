"""
Shifted Legendre basis and Galerkin infrastructure for the CL nonlinear solver.

The Galerkin solver (H&P 2017 §6) expands u and ψ as:
    u = sum_{m=0}^{J-2} sum_{k=0}^{I} A_{m,k}(t) P_m(z) cos(k l y)
    ψ = sum_{m=0}^{J}   sum_{k=0}^{I} B_{m,k}(t) P_m(z) sin(k l y)

where P_m(z) are shifted Legendre polynomials on z ∈ [-1, 0].

With I=1 harmonic and J=13 basis functions, H&P report the numerical and
asymptotic solutions are indistinguishable over l > l_cNL.

Reference: H&P (2017) §6.
"""

import numpy as np
from numpy.polynomial.legendre import legval, legder


def shifted_legendre(m: int, z: np.ndarray) -> np.ndarray:
    """
    m-th shifted Legendre polynomial P_m on z ∈ [-1, 0].

    The standard Legendre polynomial P_m is defined on [-1, 1].
    We shift via the map t = 2z + 1, so t ∈ [-1, 1] when z ∈ [-1, 0].

    Parameters:
        m: Polynomial degree [-]
        z: Evaluation points in [-1, 0] [-]

    Returns:
        P_m(z) values [-]

    Source: H&P (2017) §6; basis choice for finite-domain spectral method.
    """
    z = np.asarray(z, dtype=float)
    t = 2.0 * z + 1.0  # map [-1,0] → [-1,1]
    coeffs = np.zeros(m + 1)
    coeffs[m] = 1.0
    return legval(t, coeffs)


def shifted_legendre_deriv(m: int, order: int, z: np.ndarray) -> np.ndarray:
    """
    k-th derivative of the m-th shifted Legendre polynomial.

    Uses the chain rule: d^k/dz^k P_m(2z+1) = 2^k d^k/dt^k P_m(t)|_{t=2z+1}

    Parameters:
        m:     Polynomial degree [-]
        order: Derivative order [-]
        z:     Evaluation points in [-1, 0] [-]

    Returns:
        d^order/dz^order P_m(z) values [-]
    """
    z = np.asarray(z, dtype=float)
    t = 2.0 * z + 1.0
    coeffs = np.zeros(m + 1)
    coeffs[m] = 1.0
    # Differentiate the coefficient array `order` times
    for _ in range(order):
        coeffs = legder(coeffs)
    return (2.0 ** order) * legval(t, coeffs)


def build_basis_matrix(J: int, z: np.ndarray) -> np.ndarray:
    """
    Matrix P[m, i] = P_m(z[i]) for m = 0..J-1, evaluated on grid z.

    Parameters:
        J:  Number of basis functions [-]
        z:  Grid points [-]

    Returns:
        P: (J, n_pts) array [-]
    """
    n = len(z)
    P = np.zeros((J, n))
    for m in range(J):
        P[m] = shifted_legendre(m, z)
    return P


def build_deriv_matrix(J: int, order: int, z: np.ndarray) -> np.ndarray:
    """
    Matrix dP[m, i] = d^order P_m / dz^order at z[i].

    Parameters:
        J:     Number of basis functions [-]
        order: Derivative order [-]
        z:     Grid points [-]

    Returns:
        dP: (J, n_pts) array [-]
    """
    n = len(z)
    dP = np.zeros((J, n))
    for m in range(J):
        dP[m] = shifted_legendre_deriv(m, order, z)
    return dP


def galerkin_inner_product(
    f: np.ndarray,
    g: np.ndarray,
    z: np.ndarray,
) -> float:
    """
    Galerkin inner product ⟨f, g⟩ = ∫₋₁⁰ f(z) g(z) dz [-].

    Parameters:
        f, g: Function values on z grid [-]
        z:    Grid points [-]

    Returns:
        Inner product value [-]
    """
    return float(np.trapz(f * g, z))


def build_mass_matrix(J: int, z: np.ndarray) -> np.ndarray:
    """
    Mass matrix M[i, j] = ∫₋₁⁰ P_i(z) P_j(z) dz.

    Shifted Legendre polynomials are orthogonal with weight 1 on [-1, 0]:
        ∫₋₁⁰ P_m P_n dz = δ_{mn} / (2n+1)    (from standard orthogonality scaled by 1/2)

    Parameters:
        J: Number of basis functions [-]
        z: Integration grid [-]

    Returns:
        M: (J, J) symmetric positive-definite matrix [-]
    """
    P = build_basis_matrix(J, z)
    M = np.zeros((J, J))
    for i in range(J):
        for j in range(i, J):
            M[i, j] = np.trapz(P[i] * P[j], z)
            M[j, i] = M[i, j]
    return M


def build_stiffness_matrix(J: int, order: int, z: np.ndarray) -> np.ndarray:
    """
    Stiffness matrix K[i, j] = ∫₋₁⁰ P_i(z) d^order P_j / dz^order dz.

    Used to project the differential operators onto the basis.

    Parameters:
        J:     Number of basis functions [-]
        order: Derivative order of the right factor [-]
        z:     Integration grid [-]

    Returns:
        K: (J, J) matrix [-]
    """
    P = build_basis_matrix(J, z)
    dP = build_deriv_matrix(J, order, z)
    K = np.zeros((J, J))
    for i in range(J):
        for j in range(J):
            K[i, j] = np.trapz(P[i] * dP[j], z)
    return K


def project_function(
    f: np.ndarray,
    J: int,
    z: np.ndarray,
    M_inv: np.ndarray = None,
) -> np.ndarray:
    """
    Project f(z) onto the first J shifted Legendre basis functions.

    Returns coefficients c such that f ≈ sum_m c[m] P_m(z).

    Parameters:
        f:     Function values on z [-]
        J:     Number of basis functions [-]
        z:     Grid points [-]
        M_inv: Pre-computed inverse mass matrix (optional, computed if None)

    Returns:
        c: (J,) coefficient array [-]
    """
    P = build_basis_matrix(J, z)
    rhs = np.array([np.trapz(P[m] * f, z) for m in range(J)])
    if M_inv is None:
        M = build_mass_matrix(J, z)
        M_inv = np.linalg.inv(M)
    return M_inv @ rhs


def reconstruct(c: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    Reconstruct f(z) = sum_m c[m] P_m(z) from coefficients.

    Parameters:
        c: Coefficient array [-]
        z: Evaluation points [-]

    Returns:
        f: Reconstructed values [-]
    """
    J = len(c)
    P = build_basis_matrix(J, z)
    return c @ P
