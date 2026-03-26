"""
Scaling-law alternatives to the CL solver for Langmuir-cell geometry.

These relationships are not intended to replace the CL solver in the main
physics path. They provide an explicit alternative candidate based on the
turbulent Langmuir number La_t and a set of order-of-magnitude empirical
references for sanity checks.
"""

from __future__ import annotations

import numpy as np


def la_dependent_geometry(La_t: float, depth: float, u_star: float) -> dict:
    """
    Estimate Langmuir-cell geometry from La_t scaling laws.

    Parameters:
        La_t:   Turbulent Langmuir number [-]
        depth:  Water depth h [m]
        u_star: Water-side friction velocity [m/s]

    Returns:
        dict with scalar fields:
            downwelling_thickness:   h × La_t^{1/2} [m]
            downwelling_velocity_max: u* × La_t^{-1/3} [m/s]
            pitch:                   La_t^{1/6} [-]
            estimated_cell_width:    2 × thickness × pitch [m]

    Assumptions:
        - The literature provides proportional scalings for thickness, velocity,
          and pitch, but not a unique width closure.
        - This implementation closes width as one cell pair spanning two
          downwelling zones scaled by the LES pitch proxy.

    Source:
        McWilliams, Sullivan & Moeng (1997) for La_t; scaling summary in
        `docs/literature_review.md` §1.3.
    """
    if La_t <= 0.0:
        raise ValueError(f"La_t must be positive, got {La_t:.6g}")
    if depth <= 0.0:
        raise ValueError(f"depth must be positive, got {depth:.6g} m")
    if u_star < 0.0:
        raise ValueError(f"u_star must be non-negative, got {u_star:.6g} m/s")

    downwelling_thickness = float(depth * La_t ** 0.5)
    downwelling_velocity_max = float(u_star * La_t ** (-1.0 / 3.0))
    pitch = float(La_t ** (1.0 / 6.0))

    # Assumption AR-017: close the spanwise width from the thickness and pitch
    # scalings rather than introducing an untracked tuned prefactor.
    estimated_cell_width = float(2.0 * downwelling_thickness * pitch)

    return {
        "downwelling_thickness": downwelling_thickness,
        "downwelling_velocity_max": downwelling_velocity_max,
        "pitch": pitch,
        "estimated_cell_width": estimated_cell_width,
    }


def empirical_spacing_wind(U10: float, depth: float) -> dict:
    """
    Literature spacing ranges used as order-of-magnitude sanity checks.

    Parameters:
        U10:   Wind speed at 10 m height [m/s]
        depth: Water depth h [m]

    Returns:
        dict with scalar fields:
            faller_caponi_min_spacing: 2 × h [m]
            faller_caponi_max_spacing: 5 × h [m]
            smith_min_spacing:         2 × h [m]
            smith_max_spacing:         3 × h [m]
            composite_mid_spacing:     mean of the two literature midpoints [m]

    Notes:
        The cited sources report aspect-ratio ranges rather than a portable
        shallow-water wind-fit. This function therefore returns depth-scaled
        reference ranges, not a calibrated predictive law.

    Source:
        Faller & Caponi (1978); Smith (1992), summarised in
        `docs/literature_review.md` §3.2.
    """
    if U10 < 0.0:
        raise ValueError(f"U10 must be non-negative, got {U10:.6g} m/s")
    if depth <= 0.0:
        raise ValueError(f"depth must be positive, got {depth:.6g} m")

    faller_caponi_min_spacing = float(2.0 * depth)
    faller_caponi_max_spacing = float(5.0 * depth)
    smith_min_spacing = float(2.0 * depth)
    smith_max_spacing = float(3.0 * depth)
    composite_mid_spacing = float(
        0.5
        * (
            0.5 * (faller_caponi_min_spacing + faller_caponi_max_spacing)
            + 0.5 * (smith_min_spacing + smith_max_spacing)
        )
    )

    return {
        "faller_caponi_min_spacing": faller_caponi_min_spacing,
        "faller_caponi_max_spacing": faller_caponi_max_spacing,
        "smith_min_spacing": smith_min_spacing,
        "smith_max_spacing": smith_max_spacing,
        "composite_mid_spacing": composite_mid_spacing,
    }
