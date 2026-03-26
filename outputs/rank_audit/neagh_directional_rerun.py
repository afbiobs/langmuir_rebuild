#!/usr/bin/env python3
"""
Neagh-only directional-fetch rerun using point_summary_enriched.csv.

Assumptions made explicit in this script:
- `bathymetry` in point_summary_enriched.csv is stored as negative depth below
  the water surface, so the physical model depth is `abs(bathymetry)` [m].
- `fetch_0` ... `fetch_340` are radial fetch distances indexed by compass
  bearing in 20 degree increments, defined relative to the matched point.
  Example: `fetch_0` is the over-water path from the point toward geographic
  north, so it is the relevant fetch bin for a northerly wind.
- ERA5/Open-Meteo `wind_direction_10m` is treated as the meteorological
  direction the wind is coming from, so the default fetch lookup uses the same
  bearing. Pass `--fetch-angle-mode to` to instead use the downwind bearing.
- Each Neagh observation is matched to the nearest point-summary row by
  haversine distance because the point table does not carry observation IDs.
- The rebuilt point summary retains the annotation `category` per `group_id`,
  so matched points can be distinguished as `manual`, `stream`, or `wiggle`.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
from typing import Any

import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.evaluation.comparison import (
    _era5_cache_key,
    _load_era5_cache_frame,
    _weather_history_with_observation_sample,
    default_site_proxy_registry,
    load_observations,
    observation_time_from_proxy,
    site_from_coordinates,
)
from src.forcing import compute_forcing
from src.hydro.coarsening import disruption_check
from src.prediction.candidate_cl import analyse_candidate_cl
from src.prediction.common import build_environmental_context
from src.prediction.pipeline import (
    _disruption_history_from_weather,
    _estimate_pattern_lifetime,
    _normalise_observation_time,
    _prepare_weather_frame,
    _select_history,
)


RANK_COLUMNS = [
    "U10",
    "depth",
    "fetch_at_observation_m",
    "Ra",
    "L_inst",
    "coarsened_width",
    "comparison_spacing",
    "pattern_lifetime_s",
    "u_star_water",
    "nu_T_vertical",
    "A_H",
]


def _as_float(value: Any) -> float:
    """Convert scalar values to float, returning NaN on failure."""
    if value is None:
        return float("nan")
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _first_finite(*values: Any) -> float:
    """Return the first finite numeric value, else NaN."""
    for value in values:
        number = _as_float(value)
        if math.isfinite(number):
            return number
    return float("nan")


def _field(container: Any, name: str) -> Any:
    """Read `name` from either a dict-like object or an attribute-bearing object."""
    if isinstance(container, dict):
        return container.get(name)
    return getattr(container, name, None)


def _haversine_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance between two coordinates [m]."""
    radius = 6_371_000.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = (
        math.sin(dphi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    )
    return 2.0 * radius * math.asin(math.sqrt(a))


def _nearest_fetch_bin(direction_deg: float) -> int:
    """Round a bearing to the nearest 20 degree fetch bin."""
    return int((round((direction_deg % 360.0) / 20.0) * 20) % 360)


def _fetch_column_for_direction(direction_deg: float, angle_mode: str) -> tuple[str, int]:
    """
    Convert a wind direction to the fetch-column name.

    Parameters:
        direction_deg: Meteorological wind direction [deg]
        angle_mode:    "from" to use the incoming wind bearing, "to" to use the
                       downwind bearing [-]

    Returns:
        (column_name, fetch_bin_deg)

    Note:
        The fetch bins are point-relative azimuths. `fetch_0` means the fetch
        distance from the point toward north, so a northerly wind uses
        `fetch_0`, an easterly wind uses `fetch_90`/`fetch_100`, etc.
    """
    if angle_mode not in {"from", "to"}:
        raise ValueError(f"angle_mode must be 'from' or 'to', got {angle_mode!r}")
    bearing = direction_deg % 360.0
    if angle_mode == "to":
        bearing = (bearing + 180.0) % 360.0
    fetch_bin_deg = _nearest_fetch_bin(bearing)
    return f"fetch_{fetch_bin_deg}", fetch_bin_deg


def load_neagh_observations(observation_path: Path) -> pd.DataFrame:
    """Load only Lough Neagh observations from the standard observation table."""
    observations = load_observations(observation_path)
    registry = default_site_proxy_registry()
    mask = observations.apply(
        lambda row: site_from_coordinates(
            latitude=float(row["authoritative_lat"]),
            longitude=float(row["authoritative_lng"]),
            site_registry=registry,
        ).site_id
        == "neagh",
        axis=1,
    )
    return observations.loc[mask].copy().reset_index(drop=True)


def match_points(
    observations: pd.DataFrame,
    point_summary_path: Path,
) -> pd.DataFrame:
    """
    Match each Neagh observation to the nearest enriched point-summary row.

    Parameters:
        observations:       Neagh-only observation table [mixed]
        point_summary_path: Enriched point summary with bathymetry/fetch fields [path]

    Returns:
        Observation table augmented with matched point metadata [mixed]
    """
    points = pd.read_csv(point_summary_path)
    match_rows: list[dict[str, Any]] = []
    for _, observation in observations.iterrows():
        best_distance = float("inf")
        best_point: pd.Series | None = None
        for _, point in points.iterrows():
            distance = _haversine_m(
                float(observation["authoritative_lat"]),
                float(observation["authoritative_lng"]),
                float(point["lat"]),
                float(point["lng"]),
            )
            if distance < best_distance:
                best_distance = distance
                best_point = point
        if best_point is None:
            raise RuntimeError("Failed to match observation to a point-summary row.")

        bathymetry = _as_float(best_point["bathymetry"])
        if not math.isfinite(bathymetry) or bathymetry == 0.0:
            raise ValueError(
                f"Invalid bathymetry for matched point {best_point['group_id']}: {bathymetry}"
            )
        depth = abs(bathymetry)
        match_row = observation.to_dict()
        match_row["point_group_id"] = str(best_point["group_id"])
        match_row["annotation_category"] = str(best_point.get("category", ""))
        match_row["point_lat"] = float(best_point["lat"])
        match_row["point_lng"] = float(best_point["lng"])
        match_row["point_match_distance_m"] = float(best_distance)
        match_row["point_length_m"] = _as_float(best_point.get("length_m"))
        match_row["point_bearing_deg"] = _as_float(best_point.get("bearing_deg"))
        match_row["point_mean_spacing_m"] = _as_float(best_point.get("dist_from_prev_m"))
        match_row["point_group_index_mean"] = _as_float(best_point.get("group_index"))
        match_row["point_n_annotations"] = _as_float(best_point.get("n_points"))
        match_row["bathymetry_raw_m"] = bathymetry
        match_row["depth_m"] = depth
        for bearing in range(0, 360, 20):
            match_row[f"fetch_{bearing}"] = float(best_point[f"fetch_{bearing}"])
        match_rows.append(match_row)
    return pd.DataFrame(match_rows).sort_values("observation_id").reset_index(drop=True)


def _cache_file_for_observation(
    observation_row: pd.Series,
    cache_dir: Path,
) -> Path:
    """Resolve the ERA5 cache file used for one observation."""
    image_date = pd.Timestamp(observation_row["image_date"])
    spinup_start = (image_date.date() - pd.Timedelta(days=10)).isoformat()
    spinup_end = image_date.date().isoformat()
    cache_key = _era5_cache_key(
        float(observation_row["authoritative_lat"]),
        float(observation_row["authoritative_lng"]),
        spinup_start,
        spinup_end,
    )
    return cache_dir / f"{cache_key}.json"


def _forcing_history_with_directional_fetch(
    weather_data: pd.DataFrame,
    observation_time_utc,
    depth_m: float,
    drag_method: str,
    drift_method: str,
    angle_mode: str,
    observation_row: pd.Series,
    lookback_hours: float,
) -> tuple[list[Any], list[dict[str, Any]], pd.DataFrame]:
    """
    Build forcing history using a per-time-step fetch lookup from wind direction.

    Returns:
        (forcing_history, disruption_history, history_frame)
    """
    observation_time_utc = _normalise_observation_time(observation_time_utc)
    frame, time_col, wind_col = _prepare_weather_frame(weather_data)
    history_frame = _select_history(
        frame=frame,
        time_col=time_col,
        observation_time=observation_time_utc,
        lookback_hours=lookback_hours,
    ).copy()

    fetch_values = []
    fetch_bins = []
    for _, row in history_frame.iterrows():
        fetch_column, fetch_bin_deg = _fetch_column_for_direction(
            direction_deg=float(row["wind_direction_deg"]),
            angle_mode=angle_mode,
        )
        fetch_values.append(float(observation_row[fetch_column]))
        fetch_bins.append(fetch_bin_deg)
    history_frame["fetch_m"] = fetch_values
    history_frame["fetch_bin_deg"] = fetch_bins

    forcing_history = []
    for _, row in history_frame.iterrows():
        forcing_history.append(
            compute_forcing(
                U10=float(row[wind_col]),
                depth=depth_m,
                fetch=float(row["fetch_m"]),
                timestamp=row[time_col].to_pydatetime(),
                drag_method=drag_method,
                drift_method=drift_method,
            )
        )

    disruption_history = _disruption_history_from_weather(
        weather_history=history_frame,
        time_col=time_col,
        wind_col=wind_col,
    )
    return forcing_history, disruption_history, history_frame


def rerun_observation(
    observation_row: pd.Series,
    cache_dir: Path,
    lookback_hours: float,
    angle_mode: str,
    drag_method: str,
    drift_method: str,
) -> dict[str, Any]:
    """Rerun the CL analysis for one Neagh observation with directional fetch."""
    registry = default_site_proxy_registry()
    neagh_spec = registry["neagh"]
    observation_time_utc = observation_time_from_proxy(
        image_date=pd.Timestamp(observation_row["image_date"]),
        site_spec=neagh_spec,
    )
    cache_file = _cache_file_for_observation(observation_row, cache_dir)
    if not cache_file.exists():
        raise FileNotFoundError(
            f"ERA5 cache file missing for {observation_row['observation_id']}: {cache_file}"
        )

    era5_frame = _load_era5_cache_frame(cache_file)
    weather_data = _weather_history_with_observation_sample(
        era5_frame=era5_frame,
        observation_time_utc=observation_time_utc,
        lookback_hours=lookback_hours,
    )
    forcing_history, disruption_history, history_frame = _forcing_history_with_directional_fetch(
        weather_data=weather_data,
        observation_time_utc=observation_time_utc,
        depth_m=float(observation_row["depth_m"]),
        drag_method=drag_method,
        drift_method=drift_method,
        angle_mode=angle_mode,
        observation_row=observation_row,
        lookback_hours=lookback_hours,
    )
    pattern_lifetime = _estimate_pattern_lifetime(
        forcing_history=forcing_history,
        disruption_history=disruption_history,
        lookback_hours=lookback_hours,
    )
    disruption = disruption_check(
        forcing_history=disruption_history,
        lookback_hours=lookback_hours,
    )
    result = analyse_candidate_cl(
        forcing=forcing_history[-1],
        pattern_lifetime=pattern_lifetime,
        environmental=build_environmental_context(),
        onset_only=False,
        visible_spacing_multiplier=1.0,
        max_visible_mergers=3,
        max_visible_aspect_ratio=300.0,
        max_cell_aspect_ratio=12.0,
    )

    coarsening = result.get("coarsening", {})
    critical_result = result.get("intermediate", {}).get("critical_result", {})
    l_c = _first_finite(_field(critical_result, "l_c"), _field(critical_result, "lcNL"))
    depth_m = float(observation_row["depth_m"])
    l_inst = (
        float((2.0 * math.pi * depth_m) / l_c)
        if math.isfinite(l_c) and l_c > 0.0
        else _first_finite(coarsening.get("initial_cell_width_m"))
    )
    current_row = history_frame.iloc[-1]
    return {
        "case_id": str(observation_row["observation_id"]),
        "image_date": str(pd.Timestamp(observation_row["image_date"]).date()),
        "site_id": "neagh",
        "observed_spacing": float(observation_row["target_spacing_m"]),
        "manual_spacing": _as_float(observation_row["manual_spacing_m"]),
        "wiggle_spacing": _as_float(observation_row["wiggle_spacing_m"]),
        "measurement_method": str(observation_row["measurement_method"]),
        "point_group_id": str(observation_row["point_group_id"]),
        "annotation_category": str(observation_row.get("annotation_category", "")),
        "point_match_distance_m": float(observation_row["point_match_distance_m"]),
        "point_length_m": _as_float(observation_row["point_length_m"]),
        "point_bearing_deg": _as_float(observation_row["point_bearing_deg"]),
        "point_mean_spacing_m": _as_float(observation_row["point_mean_spacing_m"]),
        "point_group_index_mean": _as_float(observation_row["point_group_index_mean"]),
        "point_n_annotations": _as_float(observation_row["point_n_annotations"]),
        "bathymetry_raw_m": float(observation_row["bathymetry_raw_m"]),
        "depth": depth_m,
        "U10": float(result["forcing_summary"]["U10"]),
        "wind_direction_deg": float(current_row["wind_direction_deg"]),
        "fetch_bin_deg": int(current_row["fetch_bin_deg"]),
        "fetch_at_observation_m": float(current_row["fetch_m"]),
        "fetch_min_window_m": float(history_frame["fetch_m"].min()),
        "fetch_max_window_m": float(history_frame["fetch_m"].max()),
        "Ra": float(result["Ra"]),
        "L_inst": float(l_inst),
        "coarsened_width": _first_finite(coarsening.get("raw_coarsened_width_m")),
        "comparison_spacing": _first_finite(coarsening.get("visible_spacing_m")),
        "pattern_lifetime_s": float(pattern_lifetime),
        "u_star_water": float(result["forcing_summary"]["u_star_water"]),
        "nu_T_vertical": float(result["forcing_summary"]["nu_T"]),
        "A_H": float(coarsening.get("coarsening_diffusivity_m2_s", float("nan"))),
        "n_events": int(coarsening.get("n_events", 0)),
        "regime": str(result["regime"]),
        "cache_file": str(cache_file),
    }


def _spearman_summary(frame: pd.DataFrame, target_col: str) -> pd.DataFrame:
    """Compute Spearman rank correlation against the standard audit columns."""
    rows: list[dict[str, Any]] = []
    for column in RANK_COLUMNS:
        valid = frame[[target_col, column]].dropna()
        if len(valid) < 2 or valid[target_col].nunique() < 2 or valid[column].nunique() < 2:
            rho = float("nan")
            p_value = float("nan")
        else:
            rho, p_value = spearmanr(valid[target_col], valid[column])
        rows.append(
            {
                "x": column,
                "n": int(len(valid)),
                "rho": _as_float(rho),
                "p_value": _as_float(p_value),
            }
        )
    return pd.DataFrame(rows)


def _layer_flip_message(summary: pd.DataFrame) -> str:
    """Explain whether anti-correlation starts at onset, coarsening, or comparison."""
    indexed = summary.set_index("x")
    ra_rho = _as_float(indexed.at["Ra", "rho"]) if "Ra" in indexed.index else float("nan")
    onset_rho = _as_float(indexed.at["L_inst", "rho"]) if "L_inst" in indexed.index else float("nan")
    coarsened_rho = (
        _as_float(indexed.at["coarsened_width", "rho"])
        if "coarsened_width" in indexed.index
        else float("nan")
    )
    comparison_rho = (
        _as_float(indexed.at["comparison_spacing", "rho"])
        if "comparison_spacing" in indexed.index
        else float("nan")
    )
    if math.isfinite(ra_rho) and math.isfinite(onset_rho) and ra_rho >= 0.0 and onset_rho < 0.0:
        return (
            "Layer flip: Ra stays non-negative "
            f"(rho={ra_rho:.3f}), but L_inst is already anti-correlated "
            f"(rho={onset_rho:.3f})."
        )
    if math.isfinite(onset_rho) and onset_rho < 0.0:
        return (
            "Layer flip: L_inst is already anti-correlated "
            f"(rho={onset_rho:.3f}); coarsening changes it to {coarsened_rho:.3f} "
            f"and comparison spacing to {comparison_rho:.3f}."
        )
    if math.isfinite(onset_rho) and onset_rho >= 0.0 and math.isfinite(coarsened_rho) and coarsened_rho < 0.0:
        return (
            "Layer flip: onset width stays non-negative "
            f"(rho={onset_rho:.3f}), and coarsening is the first anti-correlated layer "
            f"(rho={coarsened_rho:.3f})."
        )
    if (
        math.isfinite(coarsened_rho)
        and coarsened_rho >= 0.0
        and math.isfinite(comparison_rho)
        and comparison_rho < 0.0
    ):
        return (
            "Layer flip: onset and coarsening stay non-negative, and comparison "
            f"spacing flips sign (rho={comparison_rho:.3f})."
        )
    return "Layer flip: no negative correlation appears in the onset/coarsening/comparison chain."


def print_summary(frame: pd.DataFrame) -> None:
    """Print Neagh-only rank summaries for target and wiggle spacing."""
    print(f"Neagh cases rerun: {len(frame)}")
    print(
        "Point-match distance [m]: "
        f"median={frame['point_match_distance_m'].median():.1f}, "
        f"max={frame['point_match_distance_m'].max():.1f}"
    )
    print(
        "Depth from abs(bathymetry) [m]: "
        f"min={frame['depth'].min():.3f}, max={frame['depth'].max():.3f}"
    )
    print(
        "Observation fetch [m]: "
        f"min={frame['fetch_at_observation_m'].min():.1f}, "
        f"max={frame['fetch_at_observation_m'].max():.1f}"
    )

    target_summary = _spearman_summary(frame, "observed_spacing")
    print("\nNeagh target observed_spacing")
    print(target_summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print(_layer_flip_message(target_summary))

    wiggle_frame = frame.dropna(subset=["wiggle_spacing"]).copy()
    if not wiggle_frame.empty:
        wiggle_summary = _spearman_summary(wiggle_frame, "wiggle_spacing")
        print("\nNeagh wiggle_spacing subset")
        print(wiggle_summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
        print(_layer_flip_message(wiggle_summary))


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--observations",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "observations.csv",
        help="Base observation table containing spacing targets [path]",
    )
    parser.add_argument(
        "--point-summary",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "point_summary_enriched.csv",
        help="Neagh point summary with bathymetry and directional fetch [path]",
    )
    parser.add_argument(
        "--era5-cache-dir",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "era5_cache",
        help="Directory containing matched ERA5/Open-Meteo cache JSONs [path]",
    )
    parser.add_argument(
        "--lookback-hours",
        type=float,
        default=6.0,
        help="History window used for coarsening/disruption [h]",
    )
    parser.add_argument(
        "--fetch-angle-mode",
        choices=("from", "to"),
        default="from",
        help="Map wind direction to fetch using the incoming or downwind bearing [-]",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=REPO_ROOT / "outputs" / "rank_audit" / "neagh_directional_rank_audit.csv",
        help="Path for the Neagh-only rerun audit CSV [path]",
    )
    return parser.parse_args()


def main() -> None:
    """Run the Neagh-only directional-fetch rerun and write the audit CSV."""
    args = parse_args()
    observations = load_neagh_observations(args.observations)
    matched = match_points(observations, args.point_summary)
    rows = [
        rerun_observation(
            observation_row=row,
            cache_dir=args.era5_cache_dir,
            lookback_hours=args.lookback_hours,
            angle_mode=args.fetch_angle_mode,
            drag_method="coare35",
            drift_method="webb_fox_kemper",
        )
        for _, row in matched.iterrows()
    ]
    frame = pd.DataFrame(rows).sort_values("case_id").reset_index(drop=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(frame)} rows to {args.output_csv}")
    print_summary(frame)


if __name__ == "__main__":
    main()
