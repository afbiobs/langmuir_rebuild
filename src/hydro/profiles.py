"""
Polynomial shear and drift profile representations for the CL solver.

D'(z) and U'(z) are represented as polynomials on z ∈ [-1, 0]:
    D'(z) = sum_{m=0}^{M} a_m z^m
    U'(z) = sum_{n=0}^{N} b_n z^n

The uniform case D' = U' = 1 corresponds to coeffs=[1.0] (a_0 = 1, M = 0).

Reference: Hayes & Phillips (2017) eq. (7a,b).
"""

import numpy as np
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class PolynomialProfile:
    """
    A polynomial profile D'(z) or U'(z) on z ∈ [-1, 0].

    Fields:
        coeffs: Polynomial coefficients [a_0, a_1, ..., a_M] such that
                profile(z) = sum_m coeffs[m] * z^m  [-]
        name:   Optional label for diagnostics.

    The profile is callable: profile(z) returns the value at z.
    """
    coeffs: tuple  # immutable sequence of floats
    name: str = ""

    def __call__(self, z: np.ndarray) -> np.ndarray:
        """
        Evaluate the polynomial profile at z.

        Parameters:
            z: Depth coordinate(s) [-], z ∈ [-1, 0]

        Returns:
            profile value(s) [-]
        """
        z = np.asarray(z, dtype=float)
        result = np.zeros_like(z)
        for m, a in enumerate(self.coeffs):
            result += a * z ** m
        return result

    def integrate_once(self, z: np.ndarray, C: float = 0.0) -> np.ndarray:
        """
        Indefinite integral ∫ profile(z) dz + C.

        Parameters:
            z: Depth coordinate(s) [-]
            C: Integration constant [-]

        Returns:
            Integral value(s) [-]
        """
        z = np.asarray(z, dtype=float)
        result = np.zeros_like(z)
        for m, a in enumerate(self.coeffs):
            result += a * z ** (m + 1) / (m + 1)
        return result + C

    def definite_integral(self, z_lo: float = -1.0, z_hi: float = 0.0) -> float:
        """
        Definite integral ∫_{z_lo}^{z_hi} profile(z) dz [-].

        Parameters:
            z_lo: Lower limit [-]  (default: -1)
            z_hi: Upper limit [-]  (default: 0)

        Returns:
            Definite integral value [-]
        """
        antideriv_hi = sum(a * z_hi ** (m + 1) / (m + 1) for m, a in enumerate(self.coeffs))
        antideriv_lo = sum(a * z_lo ** (m + 1) / (m + 1) for m, a in enumerate(self.coeffs))
        return float(antideriv_hi - antideriv_lo)


def polynomial_profile(coeffs: list, name: str = "") -> PolynomialProfile:
    """
    Create a polynomial profile from a list of coefficients.

    Parameters:
        coeffs: [a_0, a_1, ..., a_M] where profile(z) = sum a_m z^m  [-]
        name:   Optional label for diagnostics.

    Returns:
        PolynomialProfile instance

    Examples:
        D' = 1       → polynomial_profile([1.0])
        D' = 1 + z   → polynomial_profile([1.0, 1.0])
        U' = 1 + 2z  → polynomial_profile([1.0, 2.0])

    Source: H&P (2017) eq. (7a,b).
    """
    return PolynomialProfile(coeffs=tuple(float(c) for c in coeffs), name=name)


def integrate_product_dz(
    f: PolynomialProfile,
    g: PolynomialProfile,
    z_lo: float = -1.0,
    z_hi: float = 0.0,
) -> float:
    """
    Compute ∫_{z_lo}^{z_hi} f(z) g(z) dz analytically [-].

    Parameters:
        f, g: PolynomialProfile instances
        z_lo: Lower integration limit [-]
        z_hi: Upper integration limit [-]

    Returns:
        Definite integral value [-]

    Used extensively in the CL solver for inner products.
    """
    # Product of two polynomials is a polynomial; integrate analytically
    result = 0.0
    for m, a in enumerate(f.coeffs):
        for n, b in enumerate(g.coeffs):
            p = m + n  # power of z
            result += a * b * (z_hi ** (p + 1) - z_lo ** (p + 1)) / (p + 1)
    return float(result)


def iterated_integral(
    integrand_coeffs: np.ndarray,
    n_integrations: int,
    z_eval: np.ndarray,
    integration_constants: np.ndarray = None,
) -> np.ndarray:
    """
    Compute the n-th iterated integral of a polynomial given by coefficients.

    For the CL solver we need ∫∫ f dz^2, ∫∫∫∫ f dz^4, etc.
    Each integration introduces one constant of integration.

    Parameters:
        integrand_coeffs: Coefficient array [a_0, a_1, ..., a_M] such that
                          f(z) = sum_m a_m z^m  [-]
        n_integrations:   Number of times to integrate [-]
        z_eval:           Points at which to evaluate result [-]
        integration_constants: Array of length n_integrations; C[0] is the
                               constant for the first integration, C[1] for
                               the second, etc. Default: all zeros.

    Returns:
        Values of the n-th iterated integral at z_eval [-]

    Notes:
        Result has degree M + n_integrations.
        Integration constants are applied as: after the k-th integration,
        add constants[k-1] * z^0 (i.e. the additive constant).
        This is the "particular integral" convention used in the H&P algorithm.
    """
    z_eval = np.asarray(z_eval, dtype=float)
    if integration_constants is None:
        integration_constants = np.zeros(n_integrations)

    # Start with the integrand polynomial coefficients
    # poly[m] = coefficient of z^m
    poly = np.array(integrand_coeffs, dtype=float)

    for k in range(n_integrations):
        # Integrate: a_m z^m → a_m/(m+1) z^(m+1), then add constant
        new_poly = np.zeros(len(poly) + 1)
        for m, a in enumerate(poly):
            new_poly[m + 1] += a / (m + 1)
        new_poly[0] += integration_constants[k]
        poly = new_poly

    # Evaluate at z_eval
    result = np.zeros_like(z_eval)
    for m, a in enumerate(poly):
        result += a * z_eval ** m
    return result
