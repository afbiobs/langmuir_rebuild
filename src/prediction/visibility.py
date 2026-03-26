"""
Visibility diagnostic for Langmuir-cell surface patterns.

This module does not modify the predicted spacing. It only estimates whether a
pattern with a given convergence velocity and lifetime is likely to be visible
from satellite imagery under the prevailing wind conditions.
"""

from __future__ import annotations

import math


_CONVERGENCE_VISIBLE_THRESHOLD = 1.0e-4   # [m/s] order-of-magnitude tracer convergence
_WIND_OBSCURING_THRESHOLD = 10.0          # [m/s] roughness begins to obscure streaks


def is_pattern_visible(
    convergence_v: float,
    U10: float,
    tracer_accumulation_time: float,
    pattern_lifetime: float,
) -> dict:
    """
    Diagnose whether a Langmuir pattern is likely satellite-visible.

    Parameters:
        convergence_v:            Surface convergence velocity [m/s]
        U10:                      Wind speed at 10 m height [m/s]
        tracer_accumulation_time: Time needed to organise a visible tracer signal [s]
        pattern_lifetime:         Time the pattern is expected to persist [s]

    Returns:
        dict with scalar fields:
            visible: bool
            confidence: float on [0, 1]
            limiting_factor: str

    Assumptions:
        - Patterns are visible when convergence is strong enough to accumulate
          tracers, surface roughness is not too intense, and the lifetime exceeds
          the accumulation timescale.
        - Confidence is the product of convergence, wind, and persistence scores.

    Source:
        WP-05b specification in `docs/clean_room_rebuild_brief.md`; visibility
        acts as a diagnostic rather than a spacing correction (AR-007).
    """
    if convergence_v < 0.0:
        raise ValueError(
            f"convergence_v must be non-negative, got {convergence_v:.6g} m/s"
        )
    if U10 < 0.0:
        raise ValueError(f"U10 must be non-negative, got {U10:.6g} m/s")
    if tracer_accumulation_time <= 0.0:
        raise ValueError(
            "tracer_accumulation_time must be positive, "
            f"got {tracer_accumulation_time:.6g} s"
        )
    if pattern_lifetime < 0.0:
        raise ValueError(
            f"pattern_lifetime must be non-negative, got {pattern_lifetime:.6g} s"
        )

    convergence_score = min(
        convergence_v / _CONVERGENCE_VISIBLE_THRESHOLD, 1.0
    )
    wind_score = min(_WIND_OBSCURING_THRESHOLD / max(U10, 1.0e-12), 1.0)
    persistence_score = min(pattern_lifetime / tracer_accumulation_time, 1.0)

    scores = {
        "convergence_too_weak": convergence_score,
        "wind_obscuring": wind_score,
        "insufficient_accumulation_time": persistence_score,
    }
    limiting_factor = min(scores, key=scores.get)
    visible = all(score >= 1.0 for score in scores.values())
    if visible:
        limiting_factor = "visible"

    confidence = float(convergence_score * wind_score * persistence_score)
    confidence = float(max(0.0, min(1.0, confidence)))

    return {
        "visible": visible,
        "confidence": confidence,
        "limiting_factor": limiting_factor,
    }
