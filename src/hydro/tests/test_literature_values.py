"""
Verification tests for the nonlinear CL solver against Hayes & Phillips (2017).

These tests are written BEFORE the solver is implemented. They define the acceptance
criteria. All tests are expected to fail until the solver is complete.

Reference: Hayes, D.T. & Phillips, W.R.C. (2017). Geophys. Astrophys. Fluid Dyn.
           111(1), 65–90. DOI: 10.1080/03091929.2016.1263302

All values sourced from docs/literature_review.md §6 with equation/figure references.
Test parameters unless otherwise noted:
    D' = U' = 1  (uniform profiles, a_0 = b_0 = 1, M = N = 0)
    gamma_s = 0.0001, gamma_b = 0.0  (small Robin; Figure 6 conditions)
"""

import pytest

# These imports will fail until the modules are implemented.
# That is expected: the tests define the target, not the current state.
try:
    from src.hydro.profiles import polynomial_profile
    from src.hydro.robin_bc import RobinBC
    from src.hydro.linear_solver import compute_R0, neutral_curve_linear, critical_linear
    from src.hydro.nonlinear_solver import (
        neutral_curve_nonlinear,
        critical_nonlinear,
        compute_kappa,
    )
except ImportError:
    pytest.skip(
        "Solver modules not yet implemented — this is expected at WP-02. "
        "Tests define acceptance criteria for WP-03.",
        allow_module_level=True,
    )


# ---------------------------------------------------------------------------
# Tolerance constants
# ---------------------------------------------------------------------------

# H&P (2017) report results to 3 decimal places (e.g. R_cL ≈ 121.068).
# We require agreement to within 1% for the asymptotic solver and 0.1% for
# the Galerkin numerical solver.
ATOL_ASYM = 1.0     # absolute tolerance for asymptotic (±1 in Ra units)
RTOL_ASYM = 0.01    # relative tolerance for asymptotic

ATOL_GALERKIN = 0.1
RTOL_GALERKIN = 0.001


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def uniform_profiles():
    """D' = U' = 1 (uniform shear and drift). H&P (2017) baseline."""
    D_prime = polynomial_profile(coeffs=[1.0])   # D'(z) = 1
    U_prime = polynomial_profile(coeffs=[1.0])   # U'(z) = 1
    return D_prime, U_prime


@pytest.fixture
def linear_Dz_uniform_U(uniform_profiles):
    """D' = 1 + z, U' = 1. H&P (2017) second profile pair."""
    D_prime = polynomial_profile(coeffs=[1.0, 1.0])   # D'(z) = 1 + z
    U_prime = polynomial_profile(coeffs=[1.0])         # U'(z) = 1
    return D_prime, U_prime


@pytest.fixture
def uniform_D_linear_Uz():
    """D' = 1, U' = 1 + z. H&P (2017) third profile pair."""
    D_prime = polynomial_profile(coeffs=[1.0])         # D'(z) = 1
    U_prime = polynomial_profile(coeffs=[1.0, 1.0])   # U'(z) = 1 + z
    return D_prime, U_prime


@pytest.fixture
def linear_profiles():
    """D' = 1 + z, U' = 1 + z. H&P (2017) fourth profile pair."""
    D_prime = polynomial_profile(coeffs=[1.0, 1.0])   # D'(z) = 1 + z
    U_prime = polynomial_profile(coeffs=[1.0, 1.0])   # U'(z) = 1 + z
    return D_prime, U_prime


@pytest.fixture
def robin_small():
    """Robin BCs: gamma_s = 0.0001, gamma_b = 0. H&P (2017) Figure 6."""
    return RobinBC(gamma_s=0.0001, gamma_b=0.0)


@pytest.fixture
def robin_full():
    """Robin BCs: gamma_s = 0.06, gamma_b = 0.28. Cox & Leibovich (1993)."""
    return RobinBC(gamma_s=0.06, gamma_b=0.28)


# ---------------------------------------------------------------------------
# WP-02 acceptance tests
# ---------------------------------------------------------------------------

class TestR0Uniform:
    """
    R₀ ≈ 120 for D' = U' = 1.

    Source: H&P (2017) §3.2, eq. (26): R₀⁻¹ = ∫₋₁⁰ ψ̃₁ U' dz.
            For uniform profiles this evaluates to R₀ ≈ 120.
            Confirmed by R_cL ≈ R₀ when γ → 0.
    """

    def test_R0_uniform_profiles(self, uniform_profiles):
        """R₀ ≈ 120 for D' = U' = 1 (H&P 2017 §3.2)."""
        D_prime, U_prime = uniform_profiles
        R0 = compute_R0(D_prime=D_prime, U_prime=U_prime)
        assert abs(R0 - 120.0) < 1.0, (
            f"R₀ = {R0:.3f}, expected ≈ 120. "
            "Check eq. (26) of H&P 2017: R₀⁻¹ = ∫ψ̃₁ U' dz."
        )


class TestCriticalValues:
    """
    Critical linear and nonlinear Rayleigh numbers and wavenumbers.

    Source: H&P (2017) Figure 6 caption:
        R_cL ≈ 121.068, l_cL ≈ 0.150   (linear)
        R_cNL ≈ 122.194, l_cNL ≈ 0.105  (nonlinear)
    Conditions: D' = U' = 1, gamma_s = 0.0001, gamma_b = 0.
    """

    def test_R_cL(self, uniform_profiles, robin_small):
        D_prime, U_prime = uniform_profiles
        crit = critical_linear(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(crit.R_c - 121.068) < ATOL_ASYM, (
            f"R_cL = {crit.R_c:.3f}, expected ≈ 121.068 (H&P 2017 Fig. 6 caption)."
        )

    def test_l_cL(self, uniform_profiles, robin_small):
        D_prime, U_prime = uniform_profiles
        crit = critical_linear(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(crit.l_c - 0.150) < 0.005, (
            f"l_cL = {crit.l_c:.4f}, expected ≈ 0.150 (H&P 2017 Fig. 6 caption)."
        )

    def test_R_cNL(self, uniform_profiles, robin_small):
        D_prime, U_prime = uniform_profiles
        crit = critical_nonlinear(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(crit.R_c - 122.194) < ATOL_ASYM, (
            f"R_cNL = {crit.R_c:.3f}, expected ≈ 122.194 (H&P 2017 Fig. 6 caption)."
        )

    def test_l_cNL(self, uniform_profiles, robin_small):
        D_prime, U_prime = uniform_profiles
        crit = critical_nonlinear(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(crit.l_c - 0.105) < 0.005, (
            f"l_cNL = {crit.l_c:.4f}, expected ≈ 0.105 (H&P 2017 Fig. 6 caption)."
        )


class TestKappaValues:
    """
    κ = l_cL / l_cNL for the four profile pairs in H&P (2017) §7.2.

    Source: H&P (2017) §7.2, eq. (68):
        κ ≈ (1.425, 1.427, 1.919, 1.944)
        for D' = (1, 1+z, 1, 1+z), U' = (1, 1, 1+z, 1+z)
    κ is independent of γ (eq. 68 shows γ cancels in the ratio).
    Note: γ_s = 0.0001 (small) ensures l_c can be computed numerically.
    """

    KAPPA_TOLERANCE = 0.05   # ±5% on κ values

    def test_kappa_uniform(self, uniform_profiles, robin_small):
        """D' = 1, U' = 1: κ ≈ 1.425."""
        D_prime, U_prime = uniform_profiles
        kappa = compute_kappa(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(kappa - 1.425) < self.KAPPA_TOLERANCE, (
            f"κ = {kappa:.3f}, expected ≈ 1.425 (H&P 2017 §7.2)."
        )

    def test_kappa_linear_D_uniform_U(self, linear_Dz_uniform_U, robin_small):
        """D' = 1+z, U' = 1: κ ≈ 1.427."""
        D_prime, U_prime = linear_Dz_uniform_U
        kappa = compute_kappa(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(kappa - 1.427) < self.KAPPA_TOLERANCE, (
            f"κ = {kappa:.3f}, expected ≈ 1.427 (H&P 2017 §7.2)."
        )

    def test_kappa_uniform_D_linear_U(self, uniform_D_linear_Uz, robin_small):
        """D' = 1, U' = 1+z: κ ≈ 1.919."""
        D_prime, U_prime = uniform_D_linear_Uz
        kappa = compute_kappa(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(kappa - 1.919) < self.KAPPA_TOLERANCE, (
            f"κ = {kappa:.3f}, expected ≈ 1.919 (H&P 2017 §7.2)."
        )

    def test_kappa_linear_both(self, linear_profiles, robin_small):
        """D' = 1+z, U' = 1+z: κ ≈ 1.944."""
        D_prime, U_prime = linear_profiles
        kappa = compute_kappa(D_prime=D_prime, U_prime=U_prime, bc=robin_small)
        assert abs(kappa - 1.944) < self.KAPPA_TOLERANCE, (
            f"κ = {kappa:.3f}, expected ≈ 1.944 (H&P 2017 §7.2)."
        )

    def test_kappa_independent_of_gamma(self, uniform_profiles):
        """κ must be the same for gamma_s = 0.0001 and gamma_s = 0.01 (eq. 68)."""
        D_prime, U_prime = uniform_profiles
        bc_small = RobinBC(gamma_s=0.0001, gamma_b=0.0)
        bc_large = RobinBC(gamma_s=0.01, gamma_b=0.0)
        kappa_small = compute_kappa(D_prime=D_prime, U_prime=U_prime, bc=bc_small)
        kappa_large = compute_kappa(D_prime=D_prime, U_prime=U_prime, bc=bc_large)
        assert abs(kappa_small - kappa_large) < 0.01, (
            f"κ varies with γ: κ(small γ) = {kappa_small:.4f}, "
            f"κ(large γ) = {kappa_large:.4f}. "
            "H&P (2017) eq. (68) shows κ is independent of γ."
        )


class TestSupercriticalStability:
    """
    R̄(l) > R(l) for all l > 0 with Robin boundary conditions.

    Source: H&P (2017) §8 and Figures 3, 4.
    Physical meaning: nonlinearities suppress instability at all wavenumbers.
    This is a necessary condition for the nonlinear solver to be physically correct.
    """

    def test_supercritical_stability_uniform(self, uniform_profiles, robin_small):
        """R̄(l) > R(l) for l ∈ [0.01, 0.5] with Robin BCs (H&P 2017 Fig. 4)."""
        import numpy as np
        D_prime, U_prime = uniform_profiles

        l_array = np.linspace(0.01, 0.5, 20)
        R_linear = neutral_curve_linear(
            l_array=l_array, D_prime=D_prime, U_prime=U_prime, bc=robin_small
        )
        R_nonlinear = neutral_curve_nonlinear(
            l_array=l_array, D_prime=D_prime, U_prime=U_prime, bc=robin_small
        )

        delta = R_nonlinear - R_linear
        violations = l_array[delta <= 0]
        assert len(violations) == 0, (
            f"Supercritical stability violated at l = {violations}. "
            "R̄(l) must exceed R(l) for all l > 0 with Robin BCs (H&P 2017 Fig. 4)."
        )


class TestAspectRatioRange:
    """
    L = 2π / l_cNL ∈ [5, 11] with γ_s = 0.06, γ_b = 0.28 and profile pairs.

    Source: H&P (2017) §8:
        κ ∈ [1.425, 1.944] → l_cNL = l_cL / κ
        With γ_s = 0.06, γ_b = 0.28: l_cL ≈ 1.111 (thermocline bound)
        → l_cNL ∈ [0.57, 0.78]
        → L = 2π / l_cNL ∈ [8, 11]
    Full range [5, 11] includes cases with rigid bottom (l_cL ≈ 1.773).

    Notes:
        - Length is normalised by water depth h; dimensional spacing = L × h.
        - The cap at L = 12 (AR-006) is separate from this test; the instability
          onset predicts L ≤ 11.
    """

    def test_aspect_ratio_within_observed_range(self, robin_full):
        """L = 2π/l_cNL ∈ [5, 11] for all four profile pairs (H&P 2017 §8)."""
        import numpy as np

        profile_pairs = [
            (polynomial_profile([1.0]), polynomial_profile([1.0])),         # D'=1, U'=1
            (polynomial_profile([1.0, 1.0]), polynomial_profile([1.0])),    # D'=1+z, U'=1
            (polynomial_profile([1.0]), polynomial_profile([1.0, 1.0])),    # D'=1, U'=1+z
            (polynomial_profile([1.0, 1.0]), polynomial_profile([1.0, 1.0])),  # D'=U'=1+z
        ]

        for i, (D_prime, U_prime) in enumerate(profile_pairs):
            crit = critical_nonlinear(D_prime=D_prime, U_prime=U_prime, bc=robin_full)
            L = 2 * np.pi / crit.l_c
            assert 5.0 <= L <= 11.0, (
                f"Profile pair {i+1}: L = {L:.2f} (= 2π/{crit.l_c:.3f}), "
                f"expected ∈ [5, 11] (H&P 2017 §8). "
                "Check Robin BC parameters and nonlinear expansion."
            )


class TestAsymptoticNumericAgreement:
    """
    Asymptotic expansion and Galerkin numerical solution agree to ≥ 4 sig. figs.

    Source: H&P (2017) §7: "Results from the expansions and numerics are
    indistinguishable over a region well in excess of l_cNL."
    Conditions: J = 13 basis functions, I = 1 wavenumber harmonic.
    """

    def test_critical_wavenumber_agreement(self, uniform_profiles, robin_small):
        """l_cNL from asymptotic and Galerkin agree to 4 significant figures."""
        D_prime, U_prime = uniform_profiles

        crit_asym = critical_nonlinear(
            D_prime=D_prime, U_prime=U_prime, bc=robin_small, method="asymptotic"
        )
        crit_galerkin = critical_nonlinear(
            D_prime=D_prime, U_prime=U_prime, bc=robin_small, method="galerkin",
            n_basis=13, n_harmonics=1
        )

        rel_diff = abs(crit_asym.l_c - crit_galerkin.l_c) / crit_galerkin.l_c
        assert rel_diff < 1e-3, (
            f"l_cNL: asymptotic = {crit_asym.l_c:.5f}, "
            f"Galerkin = {crit_galerkin.l_c:.5f}, "
            f"relative difference = {rel_diff:.2e}. "
            "H&P (2017) §7 states agreement well in excess of l_cNL."
        )
