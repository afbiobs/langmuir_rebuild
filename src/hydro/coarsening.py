"""
Discrete Langmuir-cell coarsening and disruption logic.

The model follows the Y-junction merger picture described by Thorpe (2004):
coarsening occurs through discrete merger events, approximately doubling the
cell-pair width each time, and the mature width is capped by a shallow-water
aspect-ratio limit.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Mapping
from datetime import datetime, timedelta
from typing import Any


_TWO_PI_SQUARED = (2.0 * math.pi) ** 2


def coarsening_diffusivity(
    nu_T_vertical: float,
    u_star_water: float,
    depth: float,
) -> dict[str, float | str]:
    """
    Effective horizontal diffusivity for the coarsening clock [mixed units].

    Parameters:
        nu_T_vertical: Depth-averaged vertical eddy viscosity from forcing [m²/s]
        u_star_water:  Water-side friction velocity u*_w [m/s]
        depth:         Water depth h [m]

    Returns:
        dict containing:
            diffusivity_m2_s:                       Effective lateral diffusivity [m²/s]
            vertical_nu_T_m2_s:                    Input vertical viscosity [m²/s]
            lateral_mixing_length_diffusivity_m2_s: u*_w × h reference scale [m²/s]
            anisotropy_ratio:                      diffusivity / nu_T_vertical [-]
            method:                                Closure identifier [-]

    Notes:
        The forcing layer's `nu_T` is a depth-averaged vertical momentum
        diffusivity derived from the parabolic profile. The merger schedule,
        however, acts on horizontal cell width. This closure therefore uses a
        first-pass lateral mixing-length diffusivity

            A_H = u*_w × h,

        while preserving the original vertical `nu_T` in the diagnostics.
        With λ = O(2πh), the slow-time clock
        τ ~ λ² / ((2π)² A_H) collapses to the original order-of-magnitude
        turnover estimate τ = O(h / u*_w) from the engineering brief.
    """
    if nu_T_vertical < 0.0 or not math.isfinite(nu_T_vertical):
        raise ValueError(
            "nu_T_vertical must be non-negative and finite, "
            f"got {nu_T_vertical:.6g} m²/s"
        )
    if u_star_water < 0.0 or not math.isfinite(u_star_water):
        raise ValueError(
            "u_star_water must be non-negative and finite, "
            f"got {u_star_water:.6g} m/s"
        )
    if depth <= 0.0 or not math.isfinite(depth):
        raise ValueError(f"depth must be positive and finite, got {depth:.6g} m")

    lateral_mixing_length = float(u_star_water * depth)
    diffusivity = float(max(nu_T_vertical, lateral_mixing_length))
    if nu_T_vertical > 0.0:
        anisotropy_ratio = float(diffusivity / nu_T_vertical)
    elif diffusivity > 0.0:
        anisotropy_ratio = float("inf")
    else:
        anisotropy_ratio = 1.0

    return {
        "diffusivity_m2_s": diffusivity,
        "vertical_nu_T_m2_s": float(nu_T_vertical),
        "lateral_mixing_length_diffusivity_m2_s": lateral_mixing_length,
        "anisotropy_ratio": anisotropy_ratio,
        "method": "mixing_length_lateral",
    }


def coarsening_timescale(
    cell_width: float,
    nu_T: float | None = None,
    *,
    horizontal_diffusivity: float | None = None,
    slow_time_factor: float = 1.0 / _TWO_PI_SQUARED,
) -> float:
    """
    Time for one Y-junction merger event [s].

    Parameters:
        cell_width:        Current cell-pair width λ [m]
        nu_T:              Legacy alias for horizontal diffusivity [m²/s]
        horizontal_diffusivity:
                          Effective lateral diffusivity for cell merger [m²/s]
        slow_time_factor:  Dimensionless prefactor linking T = l² t to
                           τ ~ λ² / ν_T [-]

    Returns:
        tau: Coarsening timescale [s], or +inf when nu_T <= 0.

    Source:
        Hayes & Phillips (2017), Eq. (27), introduces the slow nonlinear time
        T = l² t. With λ = 2πh / l and an effective horizontal turbulent
        diffusivity ν_T, this implies a merger/interaction timescale

            τ ~ λ² / ((2π)² ν_T).

        The prefactor is kept explicit for auditability. Thorpe (2004) still
        motivates the discrete Y-junction merger picture itself.
    """
    if cell_width <= 0.0:
        raise ValueError(f"cell_width must be positive, got {cell_width:.6g} m")
    diffusivity = horizontal_diffusivity if horizontal_diffusivity is not None else nu_T
    if diffusivity is None:
        raise ValueError("A horizontal diffusivity must be provided.")
    if diffusivity <= 0.0:
        return float("inf")
    if slow_time_factor <= 0.0:
        raise ValueError(
            f"slow_time_factor must be positive, got {slow_time_factor:.6g}"
        )
    return float(slow_time_factor * cell_width ** 2 / diffusivity)


def coarsening_schedule(
    time_available: float,
    initial_width: float,
    nu_T: float | None = None,
    *,
    horizontal_diffusivity: float | None = None,
    slow_time_factor: float = 1.0 / _TWO_PI_SQUARED,
    max_events: int = 64,
) -> dict[str, float | int]:
    """
    Sequential merger schedule for width-dependent coarsening [mixed units].

    Parameters:
        time_available:    Time since onset or last reset [s]
        initial_width:     Initial cell-pair width λ₀ [m]
        nu_T:              Legacy alias for horizontal diffusivity [m²/s]
        horizontal_diffusivity:
                           Effective lateral diffusivity for cell merger [m²/s]
        slow_time_factor:  Dimensionless prefactor in τ ~ λ² / ν_T [-]
        max_events:        Numerical safeguard on completed merger count [-]

    Returns:
        dict containing:
            n_events:             Completed merger count [-]
            tau_initial_s:        First-merger timescale [s]
            tau_last_completed_s: Most recent completed merger timescale [s]
            tau_next_s:           Timescale for the next merger [s]
            time_used_s:          Time consumed by completed mergers [s]
            remaining_time_s:     Unused organisation time [s]

    Notes:
        Each completed merger doubles the active cell width, so the next merger
        becomes slower by a factor of four. This implements the self-limiting
        width dependence implied by the slow nonlinear time T = l² t.
    """
    if time_available <= 0.0:
        tau_initial = coarsening_timescale(
            cell_width=initial_width,
            nu_T=nu_T,
            horizontal_diffusivity=horizontal_diffusivity,
            slow_time_factor=slow_time_factor,
        )
        return {
            "n_events": 0,
            "tau_initial_s": float(tau_initial),
            "tau_last_completed_s": float("nan"),
            "tau_next_s": float(tau_initial),
            "time_used_s": 0.0,
            "remaining_time_s": max(float(time_available), 0.0),
        }
    if max_events < 0:
        raise ValueError(f"max_events must be non-negative, got {max_events}")

    width = float(initial_width)
    tau_initial = coarsening_timescale(
        cell_width=width,
        nu_T=nu_T,
        horizontal_diffusivity=horizontal_diffusivity,
        slow_time_factor=slow_time_factor,
    )
    if not math.isfinite(tau_initial):
        return {
            "n_events": 0,
            "tau_initial_s": float(tau_initial),
            "tau_last_completed_s": float("nan"),
            "tau_next_s": float(tau_initial),
            "time_used_s": 0.0,
            "remaining_time_s": float(time_available),
        }

    time_used = 0.0
    n_events = 0
    tau_last_completed = float("nan")
    tau_next = tau_initial

    while n_events < max_events and time_used + tau_next <= time_available + 1.0e-12:
        time_used += tau_next
        tau_last_completed = tau_next
        n_events += 1
        width *= 2.0
        tau_next = coarsening_timescale(
            cell_width=width,
            nu_T=nu_T,
            horizontal_diffusivity=horizontal_diffusivity,
            slow_time_factor=slow_time_factor,
        )

    return {
        "n_events": int(n_events),
        "tau_initial_s": float(tau_initial),
        "tau_last_completed_s": float(tau_last_completed),
        "tau_next_s": float(tau_next),
        "time_used_s": float(time_used),
        "remaining_time_s": float(max(time_available - time_used, 0.0)),
    }


def count_coarsening_events(
    time_available: float,
    initial_width: float,
    nu_T: float | None = None,
    *,
    horizontal_diffusivity: float | None = None,
    slow_time_factor: float = 1.0 / _TWO_PI_SQUARED,
    max_events: int = 64,
) -> int:
    """
    Number of discrete merger events that fit within the available time [-].

    Parameters:
        time_available:    Time since onset or last reset [s]
        initial_width:     Initial cell-pair width λ₀ [m]
        nu_T:              Legacy alias for horizontal diffusivity [m²/s]
        horizontal_diffusivity:
                           Effective lateral diffusivity for cell merger [m²/s]
        slow_time_factor:  Dimensionless prefactor in τ ~ λ² / ν_T [-]
        max_events:        Numerical safeguard on completed merger count [-]

    Returns:
        n_events: Integer count of completed merger events [-]
    """
    return int(
        coarsening_schedule(
            time_available=time_available,
            initial_width=initial_width,
            nu_T=nu_T,
            horizontal_diffusivity=horizontal_diffusivity,
            slow_time_factor=slow_time_factor,
            max_events=max_events,
        )["n_events"]
    )


def coarsened_width(
    initial_width: float,
    n_events: int,
    depth: float,
    max_aspect_ratio: float = 12.0,
) -> float:
    """
    Width after n discrete merger events [m].

    Parameters:
        initial_width:     Initial cell-pair width [m]
        n_events:          Number of merger events [-]
        depth:             Water depth h [m]
        max_aspect_ratio:  Cap on width/depth [-]

    Returns:
        width: Coarsened width [m]

    Notes:
        The aspect-ratio cap emits a warning rather than silently clipping.

    Source:
        Thorpe (2004) for discrete doubling; Marmorino et al. (2005) for the
        shallow-water aspect-ratio cap, summarised in `docs/literature_review.md`.
    """
    if initial_width <= 0.0:
        raise ValueError(f"initial_width must be positive, got {initial_width:.6g} m")
    if depth <= 0.0:
        raise ValueError(f"depth must be positive, got {depth:.6g} m")
    if n_events < 0:
        raise ValueError(f"n_events must be non-negative, got {n_events}")
    if max_aspect_ratio <= 0.0:
        raise ValueError(
            f"max_aspect_ratio must be positive, got {max_aspect_ratio:.6g}"
        )

    width = float(initial_width * (2.0 ** n_events))
    cap = float(max_aspect_ratio * depth)
    if width > cap:
        warnings.warn(
            f"Coarsened width {width:.1f} m exceeds cap {cap:.1f} m. "
            f"Capping at aspect ratio {max_aspect_ratio:.1f}.",
            stacklevel=2,
        )
        width = cap
    return width


def _entry_value(entry: Any, *names: str) -> Any:
    """Retrieve the first matching value from a mapping or object."""
    if isinstance(entry, Mapping):
        for name in names:
            if name in entry:
                return entry[name]
    else:
        for name in names:
            if hasattr(entry, name):
                return getattr(entry, name)
    return None


def _timestamp_cutoff(forcing_history: list, lookback_hours: float) -> tuple[list, Any]:
    """Trim a forcing history to the requested lookback window when timestamps exist."""
    if not forcing_history:
        return [], None
    latest_time = _entry_value(forcing_history[-1], "timestamp", "time")
    if not isinstance(latest_time, datetime):
        return list(forcing_history), None
    cutoff = latest_time - timedelta(hours=lookback_hours)
    recent = [
        entry
        for entry in forcing_history
        if (
            isinstance(_entry_value(entry, "timestamp", "time"), datetime)
            and _entry_value(entry, "timestamp", "time") >= cutoff
        )
    ]
    return recent or [forcing_history[-1]], latest_time


def _circular_difference_deg(a: float, b: float) -> float:
    """Smallest absolute angular difference [deg]."""
    diff = (a - b + 180.0) % 360.0 - 180.0
    return abs(diff)


def _disruption_event_for_latest(forcing_history: list) -> dict[str, bool]:
    """
    Evaluate whether the latest sample triggers a new disruption event.

    Parameters:
        forcing_history: Chronological episode history ending at the sample to test

    Returns:
        Dict of per-mechanism flags for the latest sample only
    """
    if not forcing_history:
        return {
            "direction_change": False,
            "low_wind_shutdown": False,
            "rapid_speed_increase": False,
            "disrupted": False,
        }

    latest_speed_raw = _entry_value(
        forcing_history[-1],
        "U10",
        "wind_speed",
        "wind_speed_10m",
    )
    if latest_speed_raw is None:
        raise ValueError("forcing_history entries must provide U10 or equivalent wind speed.")
    latest_speed = float(latest_speed_raw)
    low_wind_shutdown = latest_speed < 2.0

    if len(forcing_history) < 2:
        return {
            "direction_change": False,
            "low_wind_shutdown": bool(low_wind_shutdown),
            "rapid_speed_increase": False,
            "disrupted": bool(low_wind_shutdown),
        }

    previous_speeds = [
        float(_entry_value(entry, "U10", "wind_speed", "wind_speed_10m"))
        for entry in forcing_history[:-1]
        if _entry_value(entry, "U10", "wind_speed", "wind_speed_10m") is not None
    ]
    if len(previous_speeds) != len(forcing_history) - 1:
        raise ValueError("forcing_history entries must provide U10 or equivalent wind speed.")

    directions = [
        _entry_value(entry, "wind_direction_deg", "wind_dir_deg", "direction_deg")
        for entry in forcing_history
    ]
    direction_change = False
    if directions[-1] is not None and any(d is not None for d in directions[:-1]):
        latest_direction = float(directions[-1])
        direction_change = any(
            _circular_difference_deg(latest_direction, float(direction)) > 45.0
            for direction in directions[:-1]
            if direction is not None
        )

    rapid_speed_increase = False
    if len(previous_speeds) >= 2:
        sorted_previous = sorted(previous_speeds)
        mid = len(sorted_previous) // 2
        if len(sorted_previous) % 2 == 0:
            recent_median = 0.5 * (sorted_previous[mid - 1] + sorted_previous[mid])
        else:
            recent_median = sorted_previous[mid]
        rapid_speed_increase = (
            latest_speed - recent_median >= 3.0
            and latest_speed >= 1.5 * recent_median
        )

    disrupted = bool(direction_change or low_wind_shutdown or rapid_speed_increase)
    return {
        "direction_change": bool(direction_change),
        "low_wind_shutdown": bool(low_wind_shutdown),
        "rapid_speed_increase": bool(rapid_speed_increase),
        "disrupted": disrupted,
    }


def disruption_check(forcing_history: list, lookback_hours: float = 3.0) -> dict:
    """
    Detect forcing changes that disrupt an existing Langmuir-cell pattern.

    Parameters:
        forcing_history: Sequence of mappings or objects with at least `U10` and
                         optionally `timestamp` and `wind_direction_deg`-like fields.
        lookback_hours:  Duration of recent history to inspect [h]

    Returns:
        dict with fields:
            direction_change:     Any direction-change reset in the window [bool]
            low_wind_shutdown:    Any low-wind shutdown in the window [bool]
            rapid_speed_increase: Any rapid-increase reset in the window [bool]
            disrupted:            True if any reset occurred in the window [bool]
            reset_time:           Timestamp of the latest reset event, or None
            reset_causes:         Mechanism names for the latest reset event [-]
            event_count:          Number of reset events in the window [-]
            current_direction_deg:Latest wind direction sample [deg] or None
            current_u10:          Latest 10-m wind sample [m/s] or NaN

    Assumptions:
        - Directional disruption is triggered by a change > 45°.
        - Low-wind shutdown is triggered when the recent wind drops below 2 m/s.
        - Rapid increase is flagged when the latest U10 is at least 50% and
          3 m/s above the recent median.

    Source:
        Thorpe (2004) and Marmorino et al. (2005), summarised in
        `docs/literature_review.md` §2.2.
    """
    if lookback_hours <= 0.0:
        raise ValueError(f"lookback_hours must be positive, got {lookback_hours:.6g}")

    recent, latest_time = _timestamp_cutoff(forcing_history, lookback_hours)
    if not recent:
        return {
            "direction_change": False,
            "low_wind_shutdown": False,
            "rapid_speed_increase": False,
            "disrupted": False,
            "reset_time": None,
            "reset_causes": [],
            "event_count": 0,
            "current_direction_deg": None,
            "current_u10": float("nan"),
        }

    direction_change = False
    low_wind_shutdown = False
    rapid_speed_increase = False
    reset_time = None
    reset_causes: list[str] = []
    event_count = 0
    episode_history: list = []

    for entry in recent:
        episode_history.append(entry)
        event = _disruption_event_for_latest(episode_history)
        if event["disrupted"]:
            direction_change = bool(direction_change or event["direction_change"])
            low_wind_shutdown = bool(
                low_wind_shutdown or event["low_wind_shutdown"]
            )
            rapid_speed_increase = bool(
                rapid_speed_increase or event["rapid_speed_increase"]
            )
            reset_time = _entry_value(entry, "timestamp", "time")
            reset_causes = [
                name
                for name in (
                    "direction_change",
                    "low_wind_shutdown",
                    "rapid_speed_increase",
                )
                if event[name]
            ]
            event_count += 1
            episode_history = [entry]

    current_direction = _entry_value(
        recent[-1], "wind_direction_deg", "wind_dir_deg", "direction_deg"
    )
    current_u10 = _entry_value(recent[-1], "U10", "wind_speed", "wind_speed_10m")
    disrupted = bool(event_count > 0)

    return {
        "direction_change": bool(direction_change),
        "low_wind_shutdown": bool(low_wind_shutdown),
        "rapid_speed_increase": bool(rapid_speed_increase),
        "disrupted": disrupted,
        "reset_time": reset_time,
        "reset_causes": reset_causes,
        "event_count": int(event_count),
        "current_direction_deg": (
            float(current_direction) if current_direction is not None else None
        ),
        "current_u10": float(current_u10) if current_u10 is not None else float("nan"),
    }
