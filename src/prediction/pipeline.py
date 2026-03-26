"""
Public entry point for WP-05 single-case analysis.

`analyse_case` converts a weather time series into forcing history, estimates
the available organisation time for Langmuir cells, then dispatches to either
the CL candidate or the La-based scaling candidate.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import pandas as pd

from src.forcing import compute_forcing
from src.hydro.coarsening import disruption_check
from src.hydro.rayleigh import classify_regime
from src.hydro.robin_bc import RobinBC
from src.prediction.candidate_cl import analyse_candidate_cl
from src.prediction.candidate_scaling import analyse_candidate_scaling
from src.prediction.common import EnvironmentalContext, build_environmental_context


def _find_column(frame: pd.DataFrame, names: tuple[str, ...]) -> str:
    """Return the first matching column name."""
    for name in names:
        if name in frame.columns:
            return name
    raise ValueError(f"Expected one of columns {names}, found {tuple(frame.columns)}.")


def _normalise_observation_time(observation_time: datetime) -> datetime:
    """Ensure the observation time is timezone-aware in UTC."""
    if observation_time.tzinfo is None:
        return observation_time.replace(tzinfo=timezone.utc)
    return observation_time.astimezone(timezone.utc)


def _prepare_weather_frame(weather_data: pd.DataFrame) -> tuple[pd.DataFrame, str, str]:
    """Sort and normalise the weather dataframe used by `analyse_case`."""
    if weather_data.empty:
        raise ValueError("weather_data must contain at least one row.")

    frame = weather_data.copy()
    time_col = _find_column(frame, ("timestamp", "time", "datetime"))
    wind_col = _find_column(frame, ("U10", "wind_speed_10m", "wind_speed"))
    frame[time_col] = pd.to_datetime(frame[time_col], utc=True)
    frame = frame.sort_values(time_col).reset_index(drop=True)
    return frame, time_col, wind_col


def _select_history(
    frame: pd.DataFrame,
    time_col: str,
    observation_time: datetime,
    lookback_hours: float,
) -> pd.DataFrame:
    """Select rows at or before the analysis time within the lookback window."""
    eligible = frame.loc[frame[time_col] <= observation_time]
    if eligible.empty:
        eligible = frame.iloc[[0]]

    current_time = eligible.iloc[-1][time_col]
    lower_bound = current_time - pd.Timedelta(hours=lookback_hours)
    recent = eligible.loc[eligible[time_col] >= lower_bound]
    if recent.empty:
        recent = eligible.iloc[[-1]]
    return recent.reset_index(drop=True)


def _forcing_history_from_weather(
    weather_history: pd.DataFrame,
    time_col: str,
    wind_col: str,
    depth: float,
    fetch: float,
    drag_method: str,
    drift_method: str,
) -> list:
    """Convert a weather history into a list of forcing states."""
    forcing_history = []
    for _, row in weather_history.iterrows():
        forcing_history.append(
            compute_forcing(
                U10=float(row[wind_col]),
                depth=depth,
                fetch=fetch,
                timestamp=row[time_col].to_pydatetime(),
                drag_method=drag_method,
                drift_method=drift_method,
            )
        )
    return forcing_history


def _disruption_history_from_weather(
    weather_history: pd.DataFrame,
    time_col: str,
    wind_col: str,
) -> list[dict[str, Any]]:
    """Build the disruption-audit history from the raw weather table."""
    history = []
    direction_col = None
    for name in ("wind_direction_deg", "wind_dir_deg", "direction_deg"):
        if name in weather_history.columns:
            direction_col = name
            break

    for _, row in weather_history.iterrows():
        entry = {
            "timestamp": row[time_col].to_pydatetime(),
            "U10": float(row[wind_col]),
        }
        if direction_col is not None:
            entry["wind_direction_deg"] = float(row[direction_col])
        history.append(entry)
    return history


def _estimate_pattern_lifetime(
    forcing_history: list,
    disruption_history: list[dict[str, Any]],
    lookback_hours: float,
) -> float:
    """
    Estimate time available for organised Langmuir structures [s].

    The lifetime resets whenever the recent forcing becomes subcritical or the
    disruption logic flags a pattern reset.
    """
    if not forcing_history:
        return 0.0
    latest = forcing_history[-1]
    if classify_regime(latest.Ra) == "subcritical":
        return 0.0

    reset_time = forcing_history[0].timestamp
    for forcing in forcing_history:
        if classify_regime(forcing.Ra) == "subcritical":
            reset_time = forcing.timestamp

    disruption = disruption_check(
        forcing_history=disruption_history,
        lookback_hours=lookback_hours,
    )
    disruption_reset_time = disruption.get("reset_time")
    if (
        isinstance(disruption_reset_time, datetime)
        and disruption_reset_time > reset_time
    ):
        reset_time = disruption_reset_time

    return max((latest.timestamp - reset_time).total_seconds(), 0.0)


def analyse_case(
    weather_data: pd.DataFrame,
    observation_time: datetime,
    depth: float,
    fetch: float,
    candidate: str = "cl",
    environmental: dict | None = None,
    lookback_hours: float = 6.0,
    bc: RobinBC | None = None,
    onset_only: bool = False,
    drag_method: str = "coare35",
    drift_method: str = "webb_fox_kemper",
    visible_spacing_multiplier: float = 1.0,
    max_visible_mergers: int = 3,
    max_visible_aspect_ratio: float = 300.0,
    max_cell_aspect_ratio: float = 12.0,
) -> dict:
    """
    Full WP-05 analysis for a single observation case.

    Parameters:
        weather_data:      Time-indexed weather table with a wind-speed column
        observation_time:  Time of the observation to analyse
        depth:             Water depth h [m]
        fetch:             Wind fetch [m]
        candidate:         "cl" or "scaling"
        environmental:     Optional environmental overrides for diagnostics
        lookback_hours:    History window used for coarsening/disruption [h]
        bc:                Optional Robin BC override for the CL candidate
        onset_only:        If True, the CL candidate bypasses coarsening and hierarchy
        drag_method:       Forcing drag parameterisation
        drift_method:      Forcing Stokes-drift parameterisation
        visible_spacing_multiplier: Common observation-scale spacing multiplier [-]
        max_visible_mergers:        Maximum visible merger levels carried [-]
        max_visible_aspect_ratio:   Observation-scale safety cap [width/depth]
        max_cell_aspect_ratio:      Mechanical cell-scale cap [width/depth]

    Returns:
        Audit-friendly dict containing spacing, diagnostics, explanation, and
        the intermediate forcing/coarsening quantities.
    """
    if candidate not in {"cl", "scaling"}:
        raise ValueError(f"Unknown candidate '{candidate}'. Use 'cl' or 'scaling'.")
    if depth <= 0.0:
        raise ValueError(f"depth must be positive, got {depth:.6g} m")
    if fetch <= 0.0:
        raise ValueError(f"fetch must be positive, got {fetch:.6g} m")
    if lookback_hours <= 0.0:
        raise ValueError(f"lookback_hours must be positive, got {lookback_hours:.6g}")

    observation_time = _normalise_observation_time(observation_time)
    frame, time_col, wind_col = _prepare_weather_frame(weather_data)
    history_frame = _select_history(
        frame=frame,
        time_col=time_col,
        observation_time=observation_time,
        lookback_hours=lookback_hours,
    )
    forcing_history = _forcing_history_from_weather(
        weather_history=history_frame,
        time_col=time_col,
        wind_col=wind_col,
        depth=depth,
        fetch=fetch,
        drag_method=drag_method,
        drift_method=drift_method,
    )
    disruption_history = _disruption_history_from_weather(
        weather_history=history_frame,
        time_col=time_col,
        wind_col=wind_col,
    )

    current_forcing = forcing_history[-1]
    environmental_context = build_environmental_context(environmental)
    pattern_lifetime = _estimate_pattern_lifetime(
        forcing_history=forcing_history,
        disruption_history=disruption_history,
        lookback_hours=lookback_hours,
    )
    disruption = disruption_check(
        forcing_history=disruption_history,
        lookback_hours=lookback_hours,
    )

    if candidate == "cl":
        result = analyse_candidate_cl(
            forcing=current_forcing,
            pattern_lifetime=pattern_lifetime,
            environmental=environmental_context,
            bc=bc,
            onset_only=onset_only,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
        )
    else:
        result = analyse_candidate_scaling(
            forcing=current_forcing,
            pattern_lifetime=pattern_lifetime,
            environmental=environmental_context,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
        )

    result["observation_time"] = observation_time
    result["input_sample_time"] = current_forcing.timestamp
    result["history_points"] = len(forcing_history)
    result["pattern_lifetime_s"] = float(pattern_lifetime)
    result["onset_only"] = bool(onset_only)
    result["disruption"] = disruption
    result["environmental_context"] = environmental_context
    result["candidate_selected"] = candidate
    result["lookback_hours"] = float(lookback_hours)
    return result
