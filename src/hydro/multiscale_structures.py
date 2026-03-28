"""
Module: multiscale_structures
Layer: hydro
Purpose: Compute coexisting coherent structure scales from forcing parameters.

In fully developed Langmuir turbulence (supercritical Ra), the CL instability
sets a large-scale full-depth roll organization, but the primary LC cell width
scales with the surface wave wavelength rather than with the CL onset wavelength.

Two distinct scales coexist:
    1. Wave-tied LC spacing:  d_LC = C_wave x lambda_p  [m]
       The dominant LC cell-pair width set by interaction between Stokes drift
       and wind-driven shear. DNS by Tsai & Lu (2023) gives C_wave ~ 0.34
       (nondimensional wavenumber l_s ~ 2.9).

    2. CL onset spacing:  L_inst = 2*pi*h / l_cNL  [m]
       The full-depth roll wavelength from the nonlinear CL instability.
       This sets the large-scale organization detected by spectral methods
       but is not the directly visible convergence-line spacing.

The wave-tied scale provides the initial condition for coarsening; the CL onset
scale provides a diagnostic reference for the large-scale roll pattern.

Assumptions: AR-031, AR-032, AR-033 -- see assumptions_register.md
Key references: Tsai & Lu (2023, JFM 969 A30), Xuan & Shen (2025, JFM 1023 A4)
"""

from __future__ import annotations

import math
from dataclasses import dataclass


# Default wave-tied LC coefficient from Tsai & Lu (2023)
# Their DNS gives d_s / lambda_wave ~ 0.34 (l_s ~ 2.9) across wave steepness values.
_C_WAVE_DEFAULT = 0.34


@dataclass(frozen=True)
class MultiscaleResult:
    """
    Multi-scale coherent structure prediction from a single forcing state.

    All spatial fields in SI units [m].

    Fields:
        wave_lc_spacing_m:   Primary LC cell-pair width = C_wave x lambda_p [m]
        cl_onset_spacing_m:  Full-depth roll wavelength from CL solver [m]
        C_wave:              Wave-tied coefficient [-] (default 0.34)
        lambda_p_m:          Peak wave wavelength used [m]
        depth_m:             Water depth H [m]
        aspect_ratio_wave:   wave_lc_spacing_m / depth_m [-]
        aspect_ratio_onset:  cl_onset_spacing_m / depth_m [-]
    """

    wave_lc_spacing_m: float
    cl_onset_spacing_m: float
    C_wave: float
    lambda_p_m: float
    depth_m: float
    aspect_ratio_wave: float
    aspect_ratio_onset: float


def wave_tied_lc_spacing(
    lambda_p: float,
    C_wave: float = _C_WAVE_DEFAULT,
) -> float:
    """
    Wave-tied LC cell-pair spacing [m].

    d_LC = C_wave x lambda_p

    Parameters:
        lambda_p: Peak surface wave wavelength [m]
        C_wave:   Proportionality coefficient [-] (default 0.34)

    Returns:
        d_LC: LC cell-pair spacing [m]

    Source:
        Tsai & Lu (2023), JFM 969 A30.
        DNS at Re_tau ~ 530 shows LC pair width d_s / lambda ~ 0.34
        (nondimensional wavenumber l_s ~ 2*pi / 0.34 ~ 2.9).
        Result is robust across wave steepness ak = 0.135 to 0.22.
    """
    if lambda_p <= 0.0 or not math.isfinite(lambda_p):
        raise ValueError(
            f"lambda_p must be positive and finite, got {lambda_p:.6g} m"
        )
    if C_wave <= 0.0 or not math.isfinite(C_wave):
        raise ValueError(
            f"C_wave must be positive and finite, got {C_wave:.6g}"
        )
    return float(C_wave * lambda_p)


def compute_multiscale_result(
    lambda_p: float,
    depth: float,
    cl_onset_spacing: float,
    C_wave: float = _C_WAVE_DEFAULT,
) -> MultiscaleResult:
    """
    Assemble both structure scales into a frozen result.

    Parameters:
        lambda_p:          Peak surface wave wavelength [m]
        depth:             Water depth H [m]
        cl_onset_spacing:  CL instability onset spacing L_inst = 2*pi*h/l_c [m]
        C_wave:            Wave-tied coefficient [-] (default 0.34)

    Returns:
        MultiscaleResult with both scales and diagnostic ratios
    """
    if depth <= 0.0 or not math.isfinite(depth):
        raise ValueError(f"depth must be positive and finite, got {depth:.6g} m")
    if cl_onset_spacing <= 0.0 or not math.isfinite(cl_onset_spacing):
        raise ValueError(
            "cl_onset_spacing must be positive and finite, "
            f"got {cl_onset_spacing:.6g} m"
        )

    d_lc = wave_tied_lc_spacing(lambda_p, C_wave)

    return MultiscaleResult(
        wave_lc_spacing_m=d_lc,
        cl_onset_spacing_m=float(cl_onset_spacing),
        C_wave=float(C_wave),
        lambda_p_m=float(lambda_p),
        depth_m=float(depth),
        aspect_ratio_wave=float(d_lc / depth),
        aspect_ratio_onset=float(cl_onset_spacing / depth),
    )
