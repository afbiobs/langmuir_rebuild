"""
Robin (mixed) boundary conditions for the CL perturbation equations.

The boundary conditions for the perturbed u and ψ fields are (H&P 2017 eq. 2–3):

At z = 0 (free surface):
    u_z + γ_s u = 0       (eq. 2b)
    ψ_zz + (γ_s/2) ψ_z = 0  (eq. 2c)
    ψ = 0                 (eq. 2d)

At z = -1 (bottom):
    -u_z + γ_b u = 0      (eq. 3b)
    -ψ_zz + γ_b ψ_z = 0  (eq. 3c)
    ψ = 0                 (eq. 3d)

The γ parameters couple perturbed flow to the extra free-surface stress.
Setting γ_s = γ_b = 0 recovers Neumann (stress-free) boundary conditions,
which fail to select a preferred nonzero wavenumber at onset.

In the small-l expansion, γ is scaled as γ = l^4 γ̃ (eq. 6), so the
Robin terms enter at O(l^4). The solver stores γ in physical units and
applies the l^4 scaling internally when constructing the neutral curve.

Reference: Cox & Leibovich (1993); Hayes & Phillips (2016, 2017) §2.2.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING
import warnings

if TYPE_CHECKING:
    from src.forcing import ForcingState


@dataclass(frozen=True)
class RobinBC:
    """
    Robin boundary condition parameters for the CL perturbation equations.

    Fields:
        gamma_s: Surface Robin parameter γ_s [-]
                 Cox & Leibovich (1993) estimate: γ_s ≈ 0.06
                 H&P (2017) Fig. 6 verification: γ_s = 0.0001 (near-Neumann)
        gamma_b: Bottom Robin parameter γ_b [-]
                 Cox & Leibovich (1993) estimate: γ_b ≈ 0.28
                 H&P (2017) Fig. 6 verification: γ_b = 0.0

    Notes:
        - gamma = gamma_s + gamma_b is the combined Robin parameter.
        - l_cL = (gamma * R0 / |R*2|)^(1/4) scales as gamma^(1/4).
        - l_cNL = l_cL / kappa, where kappa depends only on D', U' profiles.
        - The Robin terms first appear at O(l^4) in the small-l expansion.

    Source: H&P (2017) eqs. (6), (66)–(70); AR-001, AR-002 in assumptions_register.md.
    """
    gamma_s: float  # [-] surface Robin parameter
    gamma_b: float  # [-] bottom Robin parameter

    def __post_init__(self):
        if self.gamma_s < 0:
            raise ValueError(f"gamma_s must be non-negative, got {self.gamma_s}")
        if self.gamma_b < 0:
            raise ValueError(f"gamma_b must be non-negative, got {self.gamma_b}")

    @property
    def gamma(self) -> float:
        """Combined Robin parameter γ = γ_s + γ_b [-]."""
        return self.gamma_s + self.gamma_b

    @property
    def is_neumann(self) -> bool:
        """True if both parameters are zero (reduces to stress-free BCs)."""
        return self.gamma_s == 0.0 and self.gamma_b == 0.0

    def gamma_tilde(self, l: float) -> tuple:
        """
        Scaled Robin parameters γ̃ = γ / l^4 used in the small-l expansion.

        Parameters:
            l: Spanwise wavenumber [-]

        Returns:
            (gamma_s_tilde, gamma_b_tilde): Scaled parameters [-]

        Source: H&P (2017) eq. (6): (γ_s, γ_b) = l^4 (γ̃_s, γ̃_b)
        """
        if l == 0.0:
            raise ValueError("Cannot compute gamma_tilde at l=0 (singularity)")
        l4 = l ** 4
        return (self.gamma_s / l4, self.gamma_b / l4)


def derive_robin_bc_from_forcing(forcing: ForcingState) -> tuple[RobinBC, dict]:
    """
    Derive Robin boundary parameters from the resolved wave field.

    Parameters:
        forcing: Clean-room forcing state with wave diagnostics

    Returns:
        (`RobinBC`, diagnostics dict)

    Closure:
        γ_s = a k_p = (H_s / 2) k_p
            Surface perturbation-stress scale from wave steepness [-]

        γ_b = γ_s / sinh(k_p h)
            Bottom perturbation coupling from the linear-wave decay of orbital
            motion to the bed [-]
    """
    gamma_s = float(forcing.wave_steepness)
    kh = float(forcing.k_p * forcing.depth)
    if kh <= 0.0 or not math.isfinite(kh):
        raise ValueError(f"k_p * depth must be positive and finite, got {kh}")

    bottom_coupling_factor = 1.0 / math.sinh(kh)
    gamma_b = float(gamma_s * bottom_coupling_factor)

    if gamma_b > gamma_s:
        warnings.warn(
            "Bottom Robin parameter exceeds the surface parameter; the wave field "
            "is strongly bed-coupled in this shallow configuration.",
            stacklevel=2,
        )

    bc = RobinBC(gamma_s=gamma_s, gamma_b=gamma_b)
    diagnostics = {
        "source": "forcing_wave_steepness_bottom_reach",
        "gamma_s_raw": gamma_s,
        "gamma_b_raw": gamma_b,
        "gamma_total_raw": float(bc.gamma),
        "wave_steepness": gamma_s,
        "bottom_coupling_factor": float(bottom_coupling_factor),
        "kh": kh,
        "lambda_p_m": float(forcing.lambda_p),
        "depth_m": float(forcing.depth),
    }
    return bc, diagnostics
