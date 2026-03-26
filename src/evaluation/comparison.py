"""
Observation-set comparison harness for WP-05 and WP-06.

The benchmark specification expects forcing to be matched from an ERA5 cache at
`data/raw/era5_cache/`, but that cache is not present in this workspace. The
default runner therefore uses an explicit provisional case builder:

1. classify each observation into a documented site bucket
2. assign site-specific depth and fetch
3. estimate a representative overpass wind speed from a simple seasonal proxy
4. evaluate both candidates at a fixed organisation time

All proxy inputs are preserved in the output tables so the full-set comparison is
auditable and can be replaced later by a forcing-matched run without changing the
evaluation logic.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, time, timedelta, timezone
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.data.era5_cache import era5_cache_key
from src.evaluation.metrics import (
    attractor_test,
    dynamic_range,
    hit_rate_within_fraction,
    pearson_correlation,
    range_coverage,
    spacing_bias,
    spacing_mae,
    spacing_mae_ratio,
    spacing_rmse,
    spacing_rmse_ratio,
    spearman_correlation,
    tail_coverage,
)
from src.evaluation.plots import (
    write_attractor_diagnostic,
    write_dynamic_range_comparison,
    write_enhancement_index_timeseries,
    write_predicted_vs_observed_plot,
    write_spacing_vs_wind_plot,
    write_tail_coverage_comparison,
)
from src.forcing import compute_forcing
from src.prediction.baseline import (
    baseline_constant,
    baseline_depth_scaled,
    baseline_linear_wind,
)
from src.prediction.candidate_cl import analyse_candidate_cl
from src.prediction.candidate_scaling import analyse_candidate_scaling
from src.prediction.common import EnvironmentalContext, build_environmental_context
from src.prediction.pipeline import analyse_case


@dataclass(frozen=True)
class SiteProxySpecification:
    """
    Provisional site-level inputs for benchmark-case construction.

    Parameters:
        site_id:                Stable site identifier [-]
        label:                  Human-readable site label [-]
        latitude_range:         Inclusive latitude bounds [deg]
        longitude_range:        Inclusive longitude bounds [deg]
        depth_m:                Representative site depth [m]
        fetch_m:                Representative dominant fetch [m]
        mean_u10_mps:           Seasonal-mean proxy wind speed [m/s]
        seasonal_amplitude_mps: Proxy seasonal wind amplitude [m/s]
        utc_offset_hours:       Local UTC offset used for overpass time [h]
        overpass_local_hour:    Assumed overpass time in local clock hours [h]
        direction_deg:          Representative wind direction [deg]
        notes:                  Audit note for the proxy choice [-]
    """

    site_id: str
    label: str
    latitude_range: tuple[float, float]
    longitude_range: tuple[float, float]
    depth_m: float
    fetch_m: float
    mean_u10_mps: float
    seasonal_amplitude_mps: float
    utc_offset_hours: float
    overpass_local_hour: float = 10.5
    direction_deg: float = 225.0
    notes: str = ""

    def contains(self, latitude: float, longitude: float) -> bool:
        """Return True when the coordinate lies inside this site's proxy bounds."""
        return (
            self.latitude_range[0] <= latitude <= self.latitude_range[1]
            and self.longitude_range[0] <= longitude <= self.longitude_range[1]
        )


@dataclass(frozen=True)
class ProvisionalObservationCase:
    """
    Fully specified proxy inputs for one observation-level comparison case.

    Parameters:
        observation_id:         Observation identifier [-]
        site_id:                Site identifier [-]
        site_label:             Human-readable site label [-]
        observed_spacing_m:     Observed windrow spacing [m]
        measurement_method:     "manual" or "wiggle" [-]
        image_date:             Observation date [-]
        latitude:               Observation latitude [deg]
        longitude:              Observation longitude [deg]
        observation_time_utc:   Assumed overpass time [UTC]
        representative_u10_mps: Representative overpass wind speed [m/s]
        depth_m:                Representative site depth [m]
        fetch_m:                Representative site fetch [m]
        pattern_lifetime_s:     Proxy organisation time [s]
        wind_direction_deg:     Representative wind direction [deg]
        environmental:          Environmental defaults or overrides
        proxy_notes:            Audit note carried into the output tables
    """

    observation_id: str
    site_id: str
    site_label: str
    observed_spacing_m: float
    measurement_method: str
    image_date: pd.Timestamp
    latitude: float
    longitude: float
    observation_time_utc: datetime
    representative_u10_mps: float
    depth_m: float
    fetch_m: float
    pattern_lifetime_s: float
    wind_direction_deg: float
    environmental: EnvironmentalContext
    proxy_notes: str
    forcing_mode: str = "provisional_site_proxy"


@dataclass(frozen=True)
class MatchedObservationCase:
    """
    Case inputs for an ERA5 cache-matched WP-05 comparison.

    Parameters mirror `ProvisionalObservationCase`, but the representative wind
    and organisation time are derived from an hourly weather history matched to
    the observation date and location.
    """

    observation_id: str
    site_id: str
    site_label: str
    observed_spacing_m: float
    measurement_method: str
    image_date: pd.Timestamp
    latitude: float
    longitude: float
    observation_time_utc: datetime
    representative_u10_mps: float
    depth_m: float
    fetch_m: float
    pattern_lifetime_s: float
    wind_direction_deg: float
    environmental: EnvironmentalContext
    proxy_notes: str
    weather_data: pd.DataFrame
    cache_file: str
    lookback_hours: float = 6.0
    forcing_mode: str = "era5_cache_matched"


def default_site_proxy_registry() -> dict[str, SiteProxySpecification]:
    """
    Provisional site registry used when matched forcing data is unavailable.

    Notes:
        - Taihu and Lough Neagh depth/fetch values come from the benchmark spec.
        - Prairie-lake and Lake-Erie-region values are explicit provisional
          placeholders that must be superseded by site-matched metadata.
    """
    return {
        "taihu": SiteProxySpecification(
            site_id="taihu",
            label="Taihu Lake",
            latitude_range=(30.0, 33.0),
            longitude_range=(119.0, 121.0),
            depth_m=2.0,
            fetch_m=30_000.0,
            mean_u10_mps=4.8,
            seasonal_amplitude_mps=1.1,
            utc_offset_hours=8.0,
            notes=(
                "Depth/fetch from benchmark spec. Wind proxy is provisional until "
                "ERA5-matched forcing is available."
            ),
        ),
        "erie": SiteProxySpecification(
            site_id="erie",
            label="Lake Erie Region",
            latitude_range=(41.0, 43.0),
            longitude_range=(-84.5, -80.0),
            depth_m=7.0,
            fetch_m=25_000.0,
            mean_u10_mps=6.2,
            seasonal_amplitude_mps=1.5,
            utc_offset_hours=-5.0,
            notes=(
                "Depth/fetch are provisional nearshore-shallow placeholders. "
                "Replace with site-matched metadata and ERA5 forcing."
            ),
        ),
        "prairie": SiteProxySpecification(
            site_id="prairie",
            label="Prairie Lakes",
            latitude_range=(47.0, 49.0),
            longitude_range=(-100.5, -97.5),
            depth_m=3.5,
            fetch_m=12_000.0,
            mean_u10_mps=5.4,
            seasonal_amplitude_mps=1.4,
            utc_offset_hours=-6.0,
            notes=(
                "Depth/fetch are provisional shallow-lake placeholders. Replace "
                "with site-matched metadata and ERA5 forcing."
            ),
        ),
        "neagh": SiteProxySpecification(
            site_id="neagh",
            label="Lough Neagh",
            latitude_range=(53.0, 56.0),
            longitude_range=(-7.5, -5.0),
            depth_m=9.0,
            fetch_m=15_000.0,
            mean_u10_mps=5.8,
            seasonal_amplitude_mps=1.6,
            utc_offset_hours=0.0,
            notes=(
                "Depth/fetch from project README and benchmark spec. Wind proxy is "
                "provisional until ERA5-matched forcing is available."
            ),
        ),
    }


def load_observations(
    observations: str | Path | pd.DataFrame = "data/raw/observations.csv",
) -> pd.DataFrame:
    """
    Load the benchmark observation table and construct the primary target spacing.

    Parameters:
        observations: CSV path or already-loaded dataframe

    Returns:
        Dataframe with explicit `target_spacing_m` and `measurement_method` columns
    """
    if isinstance(observations, pd.DataFrame):
        frame = observations.copy()
    else:
        frame = pd.read_csv(observations)

    required = {
        "observation_id",
        "image_date",
        "authoritative_lat",
        "authoritative_lng",
        "manual_spacing_m",
        "wiggle_spacing_m",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Observation table is missing required columns: {missing}")

    frame["image_date"] = pd.to_datetime(frame["image_date"])
    frame["has_manual"] = frame["manual_spacing_m"].notna()
    frame["has_wiggle"] = frame["wiggle_spacing_m"].notna()
    frame["target_spacing_m"] = frame["manual_spacing_m"].fillna(frame["wiggle_spacing_m"])
    frame["measurement_method"] = np.where(frame["has_manual"], "manual", "wiggle")

    if frame["target_spacing_m"].isna().any():
        missing_ids = frame.loc[frame["target_spacing_m"].isna(), "observation_id"].tolist()
        raise ValueError(
            "Every observation must provide a target spacing. Missing target for "
            f"{missing_ids}."
        )

    return frame.sort_values("observation_id").reset_index(drop=True)


def representative_u10_from_proxy(
    image_date: pd.Timestamp,
    site_spec: SiteProxySpecification,
    peak_day_of_year: int = 30,
    min_u10_mps: float = 2.5,
    max_u10_mps: float = 10.5,
) -> float:
    """
    Representative wind speed used by the provisional observation-set runner [m/s].

    The proxy is a clipped seasonal cosine around the site's representative mean.
    """
    doy = int(image_date.timetuple().tm_yday)
    seasonal = math.cos(2.0 * math.pi * (doy - peak_day_of_year) / 365.25)
    u10 = site_spec.mean_u10_mps + site_spec.seasonal_amplitude_mps * seasonal
    return float(np.clip(u10, min_u10_mps, max_u10_mps))


def _era5_cache_key(
    latitude: float,
    longitude: float,
    start_date: str,
    end_date: str,
) -> str:
    """Cache-key format shared with the reusable ERA5 downloader."""
    return era5_cache_key(latitude, longitude, start_date, end_date)


def _load_era5_cache_frame(cache_file: Path) -> pd.DataFrame:
    """Load one cached ERA5/Open-Meteo hourly JSON file into a dataframe."""
    with cache_file.open(encoding="utf-8") as fh:
        payload = json.load(fh)
    hourly = payload["hourly"]
    frame = pd.DataFrame(hourly)
    frame["time"] = pd.to_datetime(frame["time"], utc=True)
    return frame.set_index("time").sort_index().ffill()


def _weather_history_with_observation_sample(
    era5_frame: pd.DataFrame,
    observation_time_utc: datetime,
    lookback_hours: float,
) -> pd.DataFrame:
    """
    Convert ERA5 hourly fields to the public weather-history schema.

    A linearly interpolated sample is inserted at the assumed overpass time so
    the matched rerun does not silently snap the observation to the previous hour.
    """
    weather = pd.DataFrame(
        {
            "U10": era5_frame["wind_speed_10m"].astype(float),
            "wind_direction_deg": era5_frame["wind_direction_10m"].astype(float),
        },
        index=pd.DatetimeIndex(era5_frame.index),
    ).sort_index()

    target_time = pd.Timestamp(observation_time_utc)
    if target_time not in weather.index:
        extended_index = weather.index.union(pd.DatetimeIndex([target_time])).sort_values()
        extended = weather.reindex(extended_index)
        u_comp = -weather["U10"] * np.sin(np.deg2rad(weather["wind_direction_deg"]))
        v_comp = -weather["U10"] * np.cos(np.deg2rad(weather["wind_direction_deg"]))
        extended["U10"] = extended["U10"].interpolate(method="time").ffill().bfill()
        extended["u_comp"] = (
            u_comp.reindex(extended_index).interpolate(method="time").ffill().bfill()
        )
        extended["v_comp"] = (
            v_comp.reindex(extended_index).interpolate(method="time").ffill().bfill()
        )
        extended["wind_direction_deg"] = (
            np.rad2deg(np.arctan2(-extended["u_comp"], -extended["v_comp"])) + 360.0
        ) % 360.0
        weather = extended.drop(columns=["u_comp", "v_comp"])

    lower_bound = target_time - pd.Timedelta(hours=lookback_hours + 1.0)
    trimmed = weather.loc[(weather.index >= lower_bound) & (weather.index <= target_time)]
    if trimmed.empty:
        trimmed = weather.iloc[[-1]]
    return trimmed.reset_index(names="timestamp")


def build_matched_case(
    observation_row: pd.Series,
    cache_dir: str | Path,
    site_registry: dict[str, SiteProxySpecification] | None = None,
    lookback_hours: float = 6.0,
) -> MatchedObservationCase:
    """
    Convert an observation row into a forcing-matched comparison case.

    The weather history comes from the cached 10-day ERA5 spinup window ending on
    `image_date`, with an interpolated sample at the assumed local overpass time.
    """
    if lookback_hours <= 0.0:
        raise ValueError(f"lookback_hours must be positive, got {lookback_hours:.6g}")

    spec = site_from_coordinates(
        latitude=float(observation_row["authoritative_lat"]),
        longitude=float(observation_row["authoritative_lng"]),
        site_registry=site_registry,
    )
    image_date = pd.Timestamp(observation_row["image_date"])
    observation_time_utc = observation_time_from_proxy(image_date=image_date, site_spec=spec)
    spinup_start = (image_date.date() - timedelta(days=10)).isoformat()
    spinup_end = image_date.date().isoformat()
    cache_file = (
        Path(cache_dir)
        / (_era5_cache_key(float(observation_row["authoritative_lat"]), float(observation_row["authoritative_lng"]), spinup_start, spinup_end) + ".json")
    )
    if not cache_file.exists():
        raise FileNotFoundError(
            f"Matched ERA5 cache file not found for {observation_row['observation_id']}: "
            f"{cache_file}"
        )

    era5_frame = _load_era5_cache_frame(cache_file)
    weather_data = _weather_history_with_observation_sample(
        era5_frame=era5_frame,
        observation_time_utc=observation_time_utc,
        lookback_hours=lookback_hours,
    )
    sample_row = weather_data.loc[
        weather_data["timestamp"] == pd.Timestamp(observation_time_utc)
    ].iloc[-1]
    return MatchedObservationCase(
        observation_id=str(observation_row["observation_id"]),
        site_id=spec.site_id,
        site_label=spec.label,
        observed_spacing_m=float(observation_row["target_spacing_m"]),
        measurement_method=str(observation_row["measurement_method"]),
        image_date=image_date,
        latitude=float(observation_row["authoritative_lat"]),
        longitude=float(observation_row["authoritative_lng"]),
        observation_time_utc=observation_time_utc,
        representative_u10_mps=float(sample_row["U10"]),
        depth_m=float(spec.depth_m),
        fetch_m=float(spec.fetch_m),
        pattern_lifetime_s=float("nan"),
        wind_direction_deg=float(sample_row["wind_direction_deg"]),
        environmental=build_environmental_context(),
        proxy_notes=(
            "Hourly ERA5/Open-Meteo cache matched to observation date and "
            "location. Overpass assumed 10:30 local with an interpolated sample "
            "at the overpass time."
        ),
        weather_data=weather_data,
        cache_file=str(cache_file),
        lookback_hours=float(lookback_hours),
    )


def observation_time_from_proxy(
    image_date: pd.Timestamp,
    site_spec: SiteProxySpecification,
) -> datetime:
    """
    Assumed local overpass time converted to UTC.

    The benchmark spec notes typical morning overpass times near 10:30 local time.
    """
    hour = int(site_spec.overpass_local_hour)
    minute = int(round((site_spec.overpass_local_hour - hour) * 60.0))
    local_zone = timezone(timedelta(hours=site_spec.utc_offset_hours))
    local_dt = datetime.combine(
        image_date.date(),
        time(hour=hour, minute=minute),
        tzinfo=local_zone,
    )
    return local_dt.astimezone(timezone.utc)


def site_from_coordinates(
    latitude: float,
    longitude: float,
    site_registry: dict[str, SiteProxySpecification] | None = None,
) -> SiteProxySpecification:
    """Resolve an observation coordinate to its provisional site specification."""
    registry = site_registry or default_site_proxy_registry()
    for spec in registry.values():
        if spec.contains(latitude, longitude):
            return spec
    raise ValueError(
        "No provisional site specification matches "
        f"latitude={latitude:.5f}, longitude={longitude:.5f}."
    )


def build_provisional_case(
    observation_row: pd.Series,
    site_registry: dict[str, SiteProxySpecification] | None = None,
    pattern_lifetime_s: float = 1800.0,
) -> ProvisionalObservationCase:
    """
    Convert an observation-row record into a provisional comparison case.

    Parameters:
        observation_row:     Single observation-row series
        site_registry:       Optional proxy site registry override
        pattern_lifetime_s:  Fixed organisation time used by the proxy runner [s]

    Returns:
        ProvisionalObservationCase with explicit proxy inputs
    """
    if pattern_lifetime_s < 0.0:
        raise ValueError("pattern_lifetime_s must be non-negative.")

    spec = site_from_coordinates(
        latitude=float(observation_row["authoritative_lat"]),
        longitude=float(observation_row["authoritative_lng"]),
        site_registry=site_registry,
    )
    environmental = build_environmental_context()
    return ProvisionalObservationCase(
        observation_id=str(observation_row["observation_id"]),
        site_id=spec.site_id,
        site_label=spec.label,
        observed_spacing_m=float(observation_row["target_spacing_m"]),
        measurement_method=str(observation_row["measurement_method"]),
        image_date=pd.Timestamp(observation_row["image_date"]),
        latitude=float(observation_row["authoritative_lat"]),
        longitude=float(observation_row["authoritative_lng"]),
        observation_time_utc=observation_time_from_proxy(
            image_date=pd.Timestamp(observation_row["image_date"]),
            site_spec=spec,
        ),
        representative_u10_mps=representative_u10_from_proxy(
            image_date=pd.Timestamp(observation_row["image_date"]),
            site_spec=spec,
        ),
        depth_m=float(spec.depth_m),
        fetch_m=float(spec.fetch_m),
        pattern_lifetime_s=float(pattern_lifetime_s),
        wind_direction_deg=float(spec.direction_deg),
        environmental=environmental,
        proxy_notes=spec.notes,
    )


def build_benchmark_subsets(case_table: pd.DataFrame) -> dict[str, pd.Index]:
    """
    Construct the benchmark subsets used in the WP-05 comparison.

    Parameters:
        case_table: One row per observation case

    Returns:
        Dict mapping subset id to an observation-id index
    """
    return {
        "BM-A_full": case_table.index,
        "BM-B_cross_method": case_table.index[
            case_table["has_manual"] & case_table["has_wiggle"]
        ],
        "BM-C_narrow_spacing": case_table.index[case_table["observed_spacing_m"] < 75.0],
        "BM-D_wide_spacing": case_table.index[case_table["observed_spacing_m"] > 150.0],
        "BM-E_taihu": case_table.index[case_table["site_id"] == "taihu"],
        "BM-E_north_america": case_table.index[
            case_table["site_id"].isin(["erie", "prairie"])
        ],
        "BM-E_neagh": case_table.index[case_table["site_id"] == "neagh"],
        "proxy_low_wind": case_table.index[case_table["representative_u10_mps"] < 4.0],
        "proxy_high_wind": case_table.index[case_table["representative_u10_mps"] > 8.0],
    }


def _candidate_prediction_row(
    case: ProvisionalObservationCase | MatchedObservationCase,
    raw_observation_row: pd.Series,
    model: str,
    *,
    onset_only: bool = False,
    visible_spacing_multiplier: float = 1.0,
    max_visible_mergers: int = 3,
    max_visible_aspect_ratio: float = 300.0,
    max_cell_aspect_ratio: float = 12.0,
    drag_method: str = "coare35",
    drift_method: str = "webb_fox_kemper",
    result: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run one candidate on one comparison case and flatten the result."""
    if result is None:
        result = _run_candidate_model(
            case,
            model,
            onset_only=onset_only,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
            drag_method=drag_method,
            drift_method=drift_method,
        )

    if model not in {"cl", "scaling"}:
        raise ValueError(f"Unknown candidate model '{model}'.")

    raw_predicted_spacing = float(result["predicted_spacing_m"])
    comparison_spacing = (
        raw_predicted_spacing if math.isfinite(raw_predicted_spacing) else 0.0
    )
    enhancement = result["lc_enhancement"]
    visibility = result["visibility"]
    intermediate = result.get("intermediate", {})
    robin_closure = intermediate.get("robin_closure", {})
    pattern_lifetime = float(result.get("pattern_lifetime_s", case.pattern_lifetime_s))
    disruption = result.get("disruption", {}) if isinstance(result.get("disruption"), dict) else {}
    representative_u10 = float(result["forcing_summary"]["U10"])
    current_direction = disruption.get("current_direction_deg")
    if current_direction is None or (
        isinstance(current_direction, float) and not math.isfinite(current_direction)
    ):
        wind_direction = float(case.wind_direction_deg)
    else:
        wind_direction = float(current_direction)
    reset_causes = disruption.get("reset_causes", [])
    if isinstance(reset_causes, list):
        reset_causes_text = "|".join(str(item) for item in reset_causes)
    else:
        reset_causes_text = str(reset_causes)
    return {
        "observation_id": case.observation_id,
        "model": model,
        "model_type": "candidate",
        "site_id": case.site_id,
        "site_label": case.site_label,
        "image_date": case.image_date,
        "measurement_method": case.measurement_method,
        "observed_spacing_m": case.observed_spacing_m,
        "manual_spacing_m": float(raw_observation_row["manual_spacing_m"])
        if pd.notna(raw_observation_row["manual_spacing_m"])
        else float("nan"),
        "wiggle_spacing_m": float(raw_observation_row["wiggle_spacing_m"])
        if pd.notna(raw_observation_row["wiggle_spacing_m"])
        else float("nan"),
        "representative_u10_mps": representative_u10,
        "depth_m": case.depth_m,
        "fetch_m": case.fetch_m,
        "pattern_lifetime_s": pattern_lifetime,
        "pattern_reset_time": disruption.get("reset_time", float("nan")),
        "pattern_reset_event_count": int(disruption.get("event_count", 0)),
        "pattern_reset_causes": reset_causes_text,
        "direction_change_in_window": bool(disruption.get("direction_change", False)),
        "low_wind_shutdown_in_window": bool(
            disruption.get("low_wind_shutdown", False)
        ),
        "rapid_speed_increase_in_window": bool(
            disruption.get("rapid_speed_increase", False)
        ),
        "wind_direction_deg": wind_direction,
        "raw_predicted_spacing_m": raw_predicted_spacing,
        "comparison_spacing_m": float(comparison_spacing),
        "has_lc": bool(result["has_lc"]),
        "regime": result["regime"],
        "Ra": float(result["Ra"]),
        "La_t": float(result["La_t"]),
        "forcing_nu_T_m2_s": float(result["coarsening"]["forcing_nu_T_m2_s"]),
        "coarsening_diffusivity_m2_s": float(
            result["coarsening"]["coarsening_diffusivity_m2_s"]
        ),
        "lateral_mixing_length_diffusivity_m2_s": float(
            result["coarsening"]["lateral_mixing_length_diffusivity_m2_s"]
        ),
        "coarsening_anisotropy_ratio": float(
            result["coarsening"]["coarsening_anisotropy_ratio"]
        ),
        "coarsening_diffusivity_method": str(
            result["coarsening"]["coarsening_diffusivity_method"]
        ),
        "tau_coarsening_s": float(result["coarsening"]["tau_coarsening_s"]),
        "tau_next_coarsening_s": float(result["coarsening"]["tau_next_coarsening_s"]),
        "cap_binding": bool(result["coarsening"]["cap_binding"]),
        "n_coarsening_events": int(result["coarsening"]["n_events"]),
        "n_visible_events": int(result["coarsening"]["visible_n_events"]),
        "onset_only": bool(result.get("onset_only", onset_only)),
        "initial_cell_width_m": float(result["coarsening"]["initial_cell_width_m"])
        if math.isfinite(result["coarsening"]["initial_cell_width_m"])
        else float("nan"),
        "visible_spacing_lower_m": float(result["coarsening"]["visible_spacing_lower_m"])
        if math.isfinite(result["coarsening"]["visible_spacing_lower_m"])
        else float("nan"),
        "visible_spacing_upper_m": float(result["coarsening"]["visible_spacing_upper_m"])
        if math.isfinite(result["coarsening"]["visible_spacing_upper_m"])
        else float("nan"),
        "development_index": float(enhancement["development_index"]),
        "light_ratio": float(enhancement["light"].ratio),
        "nutrient_ratio": float(enhancement["nutrients"].ratio),
        "temperature_ratio": float(enhancement["temperature"].ratio),
        "visibility_visible": bool(visibility["visible"]),
        "visibility_confidence": float(visibility["confidence"]),
        "visibility_limiting_factor": str(visibility["limiting_factor"]),
        "selected_wavenumber_method": str(
            intermediate.get("selected_wavenumber_method", "not_applicable")
        ),
        "robin_gamma_s": float(robin_closure.get("gamma_s_raw", float("nan"))),
        "robin_gamma_b": float(robin_closure.get("gamma_b_raw", float("nan"))),
        "robin_gamma_total": float(robin_closure.get("gamma_total_raw", float("nan"))),
        "robin_closure_source": str(robin_closure.get("source", "not_reported")),
        "robin_bottom_coupling_factor": float(
            robin_closure.get("bottom_coupling_factor", float("nan"))
        ),
        "input_sample_time_utc": result.get("input_sample_time", case.observation_time_utc),
        "forcing_mode": case.forcing_mode,
        "matched_cache_file": getattr(case, "cache_file", ""),
        "proxy_notes": case.proxy_notes,
        "drag_method": str(result["forcing_summary"].get("drag_method", drag_method))
        if isinstance(result.get("forcing_summary"), dict)
        else drag_method,
        "drift_method": str(result["forcing_summary"].get("drift_method", drift_method))
        if isinstance(result.get("forcing_summary"), dict)
        else drift_method,
        "explanation": result["explanation"],
    }


def _run_candidate_model(
    case: ProvisionalObservationCase | MatchedObservationCase,
    model: str,
    *,
    onset_only: bool = False,
    visible_spacing_multiplier: float = 1.0,
    max_visible_mergers: int = 3,
    max_visible_aspect_ratio: float = 300.0,
    max_cell_aspect_ratio: float = 12.0,
    drag_method: str = "coare35",
    drift_method: str = "webb_fox_kemper",
) -> dict[str, Any]:
    """Run one candidate on one comparison case and retain full diagnostics."""
    if isinstance(case, MatchedObservationCase):
        return analyse_case(
            weather_data=case.weather_data,
            observation_time=case.observation_time_utc,
            depth=case.depth_m,
            fetch=case.fetch_m,
            candidate=model,
            environmental=asdict(case.environmental),
            lookback_hours=case.lookback_hours,
            onset_only=onset_only,
            drag_method=drag_method,
            drift_method=drift_method,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
        )
    forcing = compute_forcing(
        U10=case.representative_u10_mps,
        depth=case.depth_m,
        fetch=case.fetch_m,
        timestamp=case.observation_time_utc,
        drag_method=drag_method,
        drift_method=drift_method,
    )
    if model == "cl":
        return analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=case.pattern_lifetime_s,
            environmental=case.environmental,
            onset_only=onset_only,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
        )
    if model == "scaling":
        return analyse_candidate_scaling(
            forcing=forcing,
            pattern_lifetime=case.pattern_lifetime_s,
            environmental=case.environmental,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
        )
    raise ValueError(f"Unknown candidate model '{model}'.")


def _baseline_prediction_rows(case_table: pd.DataFrame) -> tuple[list[dict[str, Any]], dict]:
    """Fit the baselines and expand them into per-observation prediction rows."""
    observed = case_table["observed_spacing_m"].to_numpy(dtype=float)
    wind = case_table["representative_u10_mps"].to_numpy(dtype=float)
    depth = case_table["depth_m"].to_numpy(dtype=float)

    baseline_fits = {
        "baseline_constant": baseline_constant(observed),
        "baseline_linear_wind": baseline_linear_wind(wind, observed),
        "baseline_depth_scaled": baseline_depth_scaled(wind, depth, observed),
    }
    rows: list[dict[str, Any]] = []
    for model, fit in baseline_fits.items():
        for observation_id, predicted_spacing in zip(case_table.index, fit["predicted"]):
            case = case_table.loc[observation_id]
            rows.append(
                {
                    "observation_id": observation_id,
                    "model": model,
                    "model_type": "baseline",
                    "site_id": case["site_id"],
                    "site_label": case["site_label"],
                    "image_date": case["image_date"],
                    "measurement_method": case["measurement_method"],
                    "observed_spacing_m": float(case["observed_spacing_m"]),
                    "manual_spacing_m": float(case["manual_spacing_m"])
                    if pd.notna(case["manual_spacing_m"])
                    else float("nan"),
                    "wiggle_spacing_m": float(case["wiggle_spacing_m"])
                    if pd.notna(case["wiggle_spacing_m"])
                    else float("nan"),
                    "representative_u10_mps": float(case["representative_u10_mps"]),
                    "depth_m": float(case["depth_m"]),
                    "fetch_m": float(case["fetch_m"]),
                    "pattern_lifetime_s": float(case["pattern_lifetime_s"]),
                    "wind_direction_deg": float(case["wind_direction_deg"]),
                    "raw_predicted_spacing_m": float(predicted_spacing),
                    "comparison_spacing_m": float(predicted_spacing),
                    "has_lc": True,
                    "regime": "baseline",
                    "Ra": float("nan"),
                    "La_t": float("nan"),
                    "cap_binding": False,
                    "n_coarsening_events": 0,
                    "n_visible_events": 0,
                    "onset_only": False,
                    "initial_cell_width_m": float("nan"),
                    "visible_spacing_lower_m": float("nan"),
                    "visible_spacing_upper_m": float("nan"),
                    "development_index": float("nan"),
                    "light_ratio": float("nan"),
                    "nutrient_ratio": float("nan"),
                    "temperature_ratio": float("nan"),
                    "visibility_visible": False,
                    "visibility_confidence": float("nan"),
                    "visibility_limiting_factor": "baseline",
                    "selected_wavenumber_method": "baseline",
                    "robin_gamma_s": float("nan"),
                    "robin_gamma_b": float("nan"),
                    "robin_gamma_total": float("nan"),
                    "robin_closure_source": "baseline",
                    "robin_bottom_coupling_factor": float("nan"),
                    "forcing_mode": "statistical_baseline",
                    "proxy_notes": "Baseline fit on the same benchmark subset inputs.",
                    "explanation": "Statistical baseline prediction.",
                }
            )
    return rows, baseline_fits


def _metrics_row(
    model: str,
    model_type: str,
    subset_id: str,
    predicted: np.ndarray,
    observed: np.ndarray,
) -> dict[str, Any]:
    """Compute the full metrics row for one model/subset pair."""
    observed_range = (float(np.min(observed)), float(np.max(observed)))
    tail = tail_coverage(predicted, observed)
    attractor = attractor_test(predicted, observed_range=observed_range)
    return {
        "subset_id": subset_id,
        "model": model,
        "model_type": model_type,
        "n_cases": int(len(observed)),
        "rmse_m": spacing_rmse(predicted, observed),
        "mae_m": spacing_mae(predicted, observed),
        "rmse_ratio_pct": spacing_rmse_ratio(predicted, observed),
        "mae_ratio_pct": spacing_mae_ratio(predicted, observed),
        "bias_m": spacing_bias(predicted, observed),
        "pearson_r": pearson_correlation(predicted, observed),
        "spearman_rho": spearman_correlation(predicted, observed),
        "hit_rate_30pct": hit_rate_within_fraction(predicted, observed),
        "tail_coverage_low": tail["low_tail_coverage"],
        "tail_coverage_high": tail["high_tail_coverage"],
        "tail_coverage_overall": tail["overall_tail_coverage"],
        "dynamic_range_p90_p10": dynamic_range(predicted),
        "range_coverage_fraction": range_coverage(predicted, observed_range),
        "dynamic_range_pass": range_coverage(predicted, observed_range) >= 0.60,
        "attractor_pass": bool(attractor["passed"]),
        "attractor_max_fraction": attractor["max_fraction"],
        "attractor_message": attractor["message"],
    }


def _candidate_summary(
    metrics_table: pd.DataFrame,
    prediction_table: pd.DataFrame,
    model: str,
) -> dict[str, Any]:
    """Summarise the full-set behaviour of a candidate model."""
    full_metrics = (
        metrics_table.loc[
            (metrics_table["subset_id"] == "BM-A_full") & (metrics_table["model"] == model)
        ]
        .iloc[0]
        .to_dict()
    )
    rows = prediction_table.loc[prediction_table["model"] == model]
    return {
        **full_metrics,
        "cap_binding_rate": float(rows["cap_binding"].mean()),
        "visible_fraction": float(rows["visibility_visible"].mean()),
        "mean_development_index": float(rows["development_index"].mean()),
        "median_development_index": float(rows["development_index"].median()),
        "mean_light_ratio": float(rows["light_ratio"].mean()),
        "mean_nutrient_ratio": float(rows["nutrient_ratio"].mean()),
        "mean_temperature_ratio": float(rows["temperature_ratio"].mean()),
    }


def _comparison_assumptions_payload() -> dict[str, Any]:
    """Return the explicit WP-06 assumption payload recorded with outputs."""
    return {
        "shallow_lake_stratification": {
            "assumption_id": "AR-030",
            "model_parameter": "S = 0",
            "status": "Uncertain",
            "statement": (
                "Operational shallow-lake onset runs use the unstratified limit "
                "S = 0."
            ),
            "wp06_treatment": (
                "Flagged uncertainty only. Stratification is not treated as a "
                "WP-06 blocker and no ad hoc correction is applied."
            ),
        }
    }


def _select_representative_case_ids(
    case_table: pd.DataFrame,
    n_cases: int = 3,
) -> list[str]:
    """Select low-, middle-, and high-spacing representative cases."""
    ordered = case_table.sort_values("observed_spacing_m").reset_index(drop=True)
    if ordered.empty:
        return []
    if len(ordered) <= n_cases:
        return [str(value) for value in ordered["observation_id"]]

    positions = np.linspace(0, len(ordered) - 1, num=n_cases, dtype=int)
    selected: list[str] = []
    for position in positions:
        observation_id = str(ordered.iloc[int(position)]["observation_id"])
        if observation_id not in selected:
            selected.append(observation_id)
    if len(selected) == n_cases:
        return selected

    for observation_id in ordered["observation_id"]:
        observation_id = str(observation_id)
        if observation_id not in selected:
            selected.append(observation_id)
        if len(selected) == n_cases:
            break
    return selected


def _assess_production_candidate(
    metrics_table: pd.DataFrame,
    *,
    forcing_mode: str,
    provisional: bool,
) -> dict[str, Any]:
    """Assess whether the WP-06 comparison identifies a promotable winner."""
    full = metrics_table.loc[metrics_table["subset_id"] == "BM-A_full"].copy()
    full = full.set_index("model")
    baselines = ["baseline_constant", "baseline_linear_wind", "baseline_depth_scaled"]
    candidate_assessments: dict[str, Any] = {}
    dominant_candidates: list[str] = []

    for model in ("cl", "scaling"):
        other_model = "scaling" if model == "cl" else "cl"
        row = full.loc[model]
        other_row = full.loc[other_model]
        baseline_rows = full.loc[baselines]
        requirements = {
            "beats_other_candidate_rmse": float(row["rmse_m"]) < float(other_row["rmse_m"]),
            "beats_other_candidate_range_coverage": float(row["range_coverage_fraction"])
            >= float(other_row["range_coverage_fraction"]),
            "beats_other_candidate_tail_coverage": float(row["tail_coverage_overall"])
            >= float(other_row["tail_coverage_overall"]),
            "beats_all_baselines_rmse": bool(
                float(row["rmse_m"]) < float(baseline_rows["rmse_m"].min())
            ),
            "beats_all_baselines_range_coverage": bool(
                float(row["range_coverage_fraction"])
                > float(baseline_rows["range_coverage_fraction"].max())
            ),
            "beats_all_baselines_tail_coverage": bool(
                float(row["tail_coverage_overall"])
                > float(baseline_rows["tail_coverage_overall"].max())
            ),
        }
        qualifies = all(requirements.values())
        if qualifies:
            dominant_candidates.append(model)
        candidate_assessments[model] = {
            "qualifies": qualifies,
            "requirements": requirements,
            "metrics": {
                "rmse_m": float(row["rmse_m"]),
                "range_coverage_fraction": float(row["range_coverage_fraction"]),
                "tail_coverage_overall": float(row["tail_coverage_overall"]),
                "dynamic_range_p90_p10": float(row["dynamic_range_p90_p10"]),
                "attractor_pass": bool(row["attractor_pass"]),
            },
        }

    leading_candidate = min(
        ("cl", "scaling"),
        key=lambda model: float(full.loc[model, "rmse_m"]),
    )
    tentative_winner = dominant_candidates[0] if len(dominant_candidates) == 1 else None
    selected_candidate = tentative_winner if (tentative_winner and not provisional) else None

    notes = [
        "AR-030 remains explicit: the shallow-lake onset path uses S = 0.",
        "Missing stratification physics is a flagged uncertainty, not a WP-06 blocker.",
    ]
    if provisional:
        notes.append(
            "This run uses provisional site-proxy forcing, so no production promotion "
            "is allowed even if one candidate leads the metrics."
        )
    else:
        notes.append(
            "This run uses ERA5 cache-matched forcing and is eligible for production "
            "promotion if a single candidate clears the WP-06 bar."
        )

    if tentative_winner is None:
        blocker = (
            "No single physics candidate beats the other candidate and all baselines "
            "on BM-A_full across RMSE, range coverage, and overall tail coverage."
        )
    elif provisional:
        blocker = (
            f"{tentative_winner} leads the metric screen, but the forcing path is "
            "still provisional."
        )
    else:
        blocker = ""

    return {
        "subset_id": "BM-A_full",
        "forcing_mode": forcing_mode,
        "provisional": bool(provisional),
        "leading_candidate": leading_candidate,
        "selected_candidate": selected_candidate,
        "promotable": selected_candidate is not None,
        "promotion_blocker": blocker,
        "candidate_assessments": candidate_assessments,
        "notes": notes,
    }


def _gate_question(
    status: str,
    answer: str,
    basis: str,
    metrics: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a compact decision-gate question record."""
    payload = {"status": status, "answer": answer, "basis": basis}
    if metrics is not None:
        payload["metrics"] = metrics
    return payload


def _decision_gate_summary(
    metrics_table: pd.DataFrame,
    prediction_table: pd.DataFrame,
    *,
    forcing_mode: str = "provisional_site_proxy",
    provisional: bool = True,
) -> dict[str, Any]:
    """Generate the WP-05 decision-gate summary."""
    full = metrics_table.loc[metrics_table["subset_id"] == "BM-A_full"].copy()
    full = full.set_index("model")
    cl = full.loc["cl"]
    scaling = full.loc["scaling"]
    baseline_best_rmse = float(
        full.loc[
            ["baseline_constant", "baseline_linear_wind", "baseline_depth_scaled"],
            "rmse_m",
        ].min()
    )

    if (
        cl["rmse_m"] < scaling["rmse_m"]
        and cl["mae_m"] <= scaling["mae_m"]
        and bool(cl["attractor_pass"]) >= bool(scaling["attractor_pass"])
    ):
        q1 = _gate_question(
            status="pass",
            answer="Yes",
            basis="CL leads the scaling candidate on fit without a worse structural check.",
            metrics={
                "cl_rmse_m": float(cl["rmse_m"]),
                "scaling_rmse_m": float(scaling["rmse_m"]),
                "cl_mae_m": float(cl["mae_m"]),
                "scaling_mae_m": float(scaling["mae_m"]),
            },
        )
    elif (
        cl["rmse_m"] >= scaling["rmse_m"]
        or cl["mae_m"] > scaling["mae_m"]
        or (not bool(cl["attractor_pass"]) and bool(scaling["attractor_pass"]))
    ):
        q1 = _gate_question(
            status="fail",
            answer="No",
            basis=(
                "CL does not clearly beat the scaling candidate once fit and "
                "structural checks are considered together."
            ),
            metrics={
                "cl_rmse_m": float(cl["rmse_m"]),
                "scaling_rmse_m": float(scaling["rmse_m"]),
                "cl_attractor_pass": bool(cl["attractor_pass"]),
                "scaling_attractor_pass": bool(scaling["attractor_pass"]),
            },
        )
    else:
        q1 = _gate_question(
            status="uncertain",
            answer="Mixed",
            basis="Fit and structural diagnostics disagree on the ranking.",
        )

    both_beat_baselines = bool(
        cl["rmse_m"] < baseline_best_rmse and scaling["rmse_m"] < baseline_best_rmse
    )
    q2 = _gate_question(
        status="pass" if both_beat_baselines else "fail",
        answer="Yes" if both_beat_baselines else "No",
        basis=(
            "Both candidates beat the best baseline on RMSE."
            if both_beat_baselines
            else "At least one baseline still fits the full dataset better than both candidates."
        ),
        metrics={
            "cl_rmse_m": float(cl["rmse_m"]),
            "scaling_rmse_m": float(scaling["rmse_m"]),
            "best_baseline_rmse_m": baseline_best_rmse,
        },
    )

    both_pass_attractor = bool(cl["attractor_pass"] and scaling["attractor_pass"])
    q3 = _gate_question(
        status="pass" if both_pass_attractor else "fail",
        answer="Yes" if both_pass_attractor else "No",
        basis=(
            "Both candidates avoid the attractor failure mode on the full dataset."
            if both_pass_attractor
            else "At least one candidate still concentrates predictions in too small a band."
        ),
        metrics={
            "cl_attractor_pass": bool(cl["attractor_pass"]),
            "scaling_attractor_pass": bool(scaling["attractor_pass"]),
            "cl_attractor_max_fraction": float(cl["attractor_max_fraction"]),
            "scaling_attractor_max_fraction": float(scaling["attractor_max_fraction"]),
        },
    )

    candidate_rows = prediction_table.loc[prediction_table["model"].isin(["cl", "scaling"])]
    ratios_positive = bool(
        (candidate_rows["light_ratio"] > 0.0).all()
        and (candidate_rows["nutrient_ratio"] > 0.0).all()
        and (candidate_rows["temperature_ratio"] > 0.0).all()
    )
    q4 = _gate_question(
        status="uncertain",
        answer="Provisionally yes",
        basis=(
            "The enhancement ratios remain positive and interpretable, but the "
            "comparison still uses default environmental inputs and simplified "
            "site metadata."
            if ratios_positive
            else "One or more enhancement ratios were non-physical."
        ),
        metrics={
            "mean_light_ratio": float(candidate_rows["light_ratio"].mean()),
            "mean_nutrient_ratio": float(candidate_rows["nutrient_ratio"].mean()),
            "mean_temperature_ratio": float(candidate_rows["temperature_ratio"].mean()),
        },
    )

    q5 = _gate_question(
        status="blocked",
        answer="Cannot assess",
        basis=(
            "The repository contains spacing observations only. There is no bloom-"
            "condition target with which to test development-index correlation."
        ),
    )

    return {
        "provisional": bool(provisional),
        "forcing_mode": forcing_mode,
        "can_advance_to_wp06": False,
        "question_1_cl_vs_scaling": q1,
        "question_2_both_beat_baselines": q2,
        "question_3_both_pass_attractor": q3,
        "question_4_enhancement_sensible": q4,
        "question_5_development_index_correlation": q5,
    }


def _json_ready(value: Any) -> Any:
    """Convert nested results to JSON-safe builtins."""
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {key: _json_ready(val) for key, val in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return str(value)


def _float_diff(left: float, right: float) -> float:
    """Return an absolute difference that treats paired NaNs as equal."""
    if math.isnan(left) and math.isnan(right):
        return 0.0
    if math.isnan(left) or math.isnan(right):
        return float("inf")
    return abs(left - right)


def _provisional_case_timeline(
    case: ProvisionalObservationCase,
    *,
    onset_only: bool,
    visible_spacing_multiplier: float,
    max_visible_mergers: int,
    max_visible_aspect_ratio: float,
    max_cell_aspect_ratio: float,
    drag_method: str,
    drift_method: str,
) -> pd.DataFrame:
    """Build a proxy enhancement timeline for a provisional comparison case."""
    forcing = compute_forcing(
        U10=case.representative_u10_mps,
        depth=case.depth_m,
        fetch=case.fetch_m,
        timestamp=case.observation_time_utc,
        drag_method=drag_method,
        drift_method=drift_method,
    )
    if case.pattern_lifetime_s <= 0.0:
        lifetime_grid = np.array([0.0])
    else:
        lifetime_grid = np.linspace(0.0, case.pattern_lifetime_s, num=6)
    start_time = case.observation_time_utc - timedelta(seconds=float(case.pattern_lifetime_s))

    rows: list[dict[str, Any]] = []
    for elapsed_s in lifetime_grid:
        timestamp = start_time + timedelta(seconds=float(elapsed_s))
        cl_result = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=float(elapsed_s),
            environmental=case.environmental,
            onset_only=onset_only,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
        )
        scaling_result = analyse_candidate_scaling(
            forcing=forcing,
            pattern_lifetime=float(elapsed_s),
            environmental=case.environmental,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
        )
        for model, result in (("cl", cl_result), ("scaling", scaling_result)):
            rows.append(
                {
                    "timestamp": timestamp,
                    "model": model,
                    "development_index": float(result["lc_enhancement"]["development_index"]),
                    "predicted_spacing_m": float(result["predicted_spacing_m"])
                    if math.isfinite(result["predicted_spacing_m"])
                    else float("nan"),
                    "timeline_mode": "proxy_lifetime_sweep",
                }
            )
    return pd.DataFrame(rows)


def _matched_case_timeline(
    case: MatchedObservationCase,
    *,
    onset_only: bool,
    visible_spacing_multiplier: float,
    max_visible_mergers: int,
    max_visible_aspect_ratio: float,
    max_cell_aspect_ratio: float,
    drag_method: str,
    drift_method: str,
) -> pd.DataFrame:
    """Build a public-pipeline enhancement timeline for a matched comparison case."""
    rows: list[dict[str, Any]] = []
    timestamps = pd.to_datetime(case.weather_data["timestamp"], utc=True)
    for timestamp in timestamps:
        for model in ("cl", "scaling"):
            result = analyse_case(
                weather_data=case.weather_data,
                observation_time=timestamp.to_pydatetime(),
                depth=case.depth_m,
                fetch=case.fetch_m,
                candidate=model,
                environmental=asdict(case.environmental),
                lookback_hours=case.lookback_hours,
                onset_only=bool(onset_only and model == "cl"),
                drag_method=drag_method,
                drift_method=drift_method,
                visible_spacing_multiplier=visible_spacing_multiplier,
                max_visible_mergers=max_visible_mergers,
                max_visible_aspect_ratio=max_visible_aspect_ratio,
                max_cell_aspect_ratio=max_cell_aspect_ratio,
            )
            rows.append(
                {
                    "timestamp": timestamp,
                    "model": model,
                    "development_index": float(result["lc_enhancement"]["development_index"]),
                    "predicted_spacing_m": float(result["predicted_spacing_m"])
                    if math.isfinite(result["predicted_spacing_m"])
                    else float("nan"),
                    "timeline_mode": "matched_weather_history",
                }
            )
    return pd.DataFrame(rows)


def _representative_case_timeline(
    case: ProvisionalObservationCase | MatchedObservationCase,
    *,
    onset_only: bool,
    visible_spacing_multiplier: float,
    max_visible_mergers: int,
    max_visible_aspect_ratio: float,
    max_cell_aspect_ratio: float,
    drag_method: str,
    drift_method: str,
) -> pd.DataFrame:
    """Build the development-index timeline for one representative case."""
    if isinstance(case, MatchedObservationCase):
        return _matched_case_timeline(
            case,
            onset_only=onset_only,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
            drag_method=drag_method,
            drift_method=drift_method,
        )
    return _provisional_case_timeline(
        case,
        onset_only=onset_only,
        visible_spacing_multiplier=visible_spacing_multiplier,
        max_visible_mergers=max_visible_mergers,
        max_visible_aspect_ratio=max_visible_aspect_ratio,
        max_cell_aspect_ratio=max_cell_aspect_ratio,
        drag_method=drag_method,
        drift_method=drift_method,
    )


def _case_diagnostic_payload(
    case_row: pd.Series,
    model_rows: pd.DataFrame,
    candidate_diagnostics: dict[str, Any],
    baseline_coefficients: dict[str, Any],
    scenario: dict[str, Any],
    *,
    forcing_mode: str,
    provisional: bool,
    representative_case_ids: list[str],
) -> dict[str, Any]:
    """Build the per-case JSON payload written into `case_diagnostics/`."""
    models: dict[str, Any] = {}
    for _, row in model_rows.iterrows():
        model = str(row["model"])
        models[model] = {"prediction_row": _json_ready(row.to_dict())}
        if model in candidate_diagnostics:
            models[model]["full_diagnostic"] = _json_ready(candidate_diagnostics[model])
        elif model in baseline_coefficients:
            models[model]["baseline_coefficients"] = _json_ready(
                baseline_coefficients[model]
            )

    return {
        "observation_id": str(case_row["observation_id"]),
        "case_input": _json_ready(case_row.to_dict()),
        "comparison_context": {
            "forcing_mode": forcing_mode,
            "provisional": bool(provisional),
            "representative_case": str(case_row["observation_id"]) in representative_case_ids,
            "scenario": _json_ready(scenario),
        },
        "assumptions": _comparison_assumptions_payload(),
        "models": models,
    }


def verify_production_candidate_reproduction(
    comparison_result: dict[str, Any],
    candidate: str,
) -> dict[str, Any]:
    """
    Re-run the public production pipeline for one candidate and compare outputs.

    Parameters:
        comparison_result: Result returned by `run_wp05_comparison` [-]
        candidate:         Candidate identifier, either "cl" or "scaling" [-]

    Returns:
        Dict with verification status and maximum metric differences [-]
    """
    if candidate not in {"cl", "scaling"}:
        raise ValueError(f"Unknown candidate '{candidate}'.")

    cases = comparison_result.get("cases")
    if not cases:
        raise ValueError("comparison_result does not contain the comparison cases.")
    if any(not isinstance(case, MatchedObservationCase) for case in cases):
        return {
            "status": "blocked",
            "candidate": candidate,
            "reason": (
                "Exact production reproduction requires matched weather histories. "
                "This comparison run was not fully matched."
            ),
        }

    prediction_rows = (
        comparison_result["prediction_table"]
        .loc[comparison_result["prediction_table"]["model"] == candidate]
        .set_index("observation_id", drop=False)
    )
    scenario = comparison_result["scenario"]
    spacing_diffs: list[float] = []
    ra_diffs: list[float] = []
    development_diffs: list[float] = []
    regime_mismatches: list[str] = []

    for case in cases:
        rerun = _run_candidate_model(
            case,
            candidate,
            onset_only=bool(scenario["onset_only"] and candidate == "cl"),
            visible_spacing_multiplier=float(scenario["visible_spacing_multiplier"]),
            max_visible_mergers=int(scenario["max_visible_mergers"]),
            max_visible_aspect_ratio=float(scenario["max_visible_aspect_ratio"]),
            max_cell_aspect_ratio=float(scenario["max_cell_aspect_ratio"]),
            drag_method=str(scenario["drag_method"]),
            drift_method=str(scenario["drift_method"]),
        )
        original = prediction_rows.loc[case.observation_id]
        spacing_diffs.append(
            _float_diff(
                float(rerun["predicted_spacing_m"]),
                float(original["raw_predicted_spacing_m"]),
            )
        )
        ra_diffs.append(_float_diff(float(rerun["Ra"]), float(original["Ra"])))
        development_diffs.append(
            _float_diff(
                float(rerun["lc_enhancement"]["development_index"]),
                float(original["development_index"]),
            )
        )
        if str(rerun["regime"]) != str(original["regime"]):
            regime_mismatches.append(str(case.observation_id))

    max_spacing_diff = max(spacing_diffs, default=0.0)
    max_ra_diff = max(ra_diffs, default=0.0)
    max_development_diff = max(development_diffs, default=0.0)
    passed = (
        max_spacing_diff <= 1.0e-9
        and max_ra_diff <= 1.0e-9
        and max_development_diff <= 1.0e-9
        and not regime_mismatches
    )
    return {
        "status": "passed" if passed else "failed",
        "candidate": candidate,
        "n_cases": len(cases),
        "max_abs_spacing_diff_m": float(max_spacing_diff),
        "max_abs_ra_diff": float(max_ra_diff),
        "max_abs_development_index_diff": float(max_development_diff),
        "regime_mismatches": regime_mismatches,
    }


def _finite_summary(values: pd.Series) -> dict[str, float]:
    """Compact min/percentile/max summary with NaNs removed."""
    finite = values.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {
            "count": 0,
            "min": float("nan"),
            "p10": float("nan"),
            "median": float("nan"),
            "p90": float("nan"),
            "max": float("nan"),
        }
    return {
        "count": int(finite.size),
        "min": float(np.min(finite)),
        "p10": float(np.quantile(finite, 0.10)),
        "median": float(np.quantile(finite, 0.50)),
        "p90": float(np.quantile(finite, 0.90)),
        "max": float(np.max(finite)),
    }


def _coarsening_closure_audit_payload(prediction_table: pd.DataFrame) -> dict[str, Any]:
    """Summarise forcing ν_T and merger-clock diffusivity across candidate runs."""
    candidate_rows = prediction_table.loc[
        prediction_table["model"].isin(["cl", "scaling"])
    ].copy()
    payload: dict[str, Any] = {
        "models": {},
        "notes": [
            "forcing_nu_T_m2_s is the depth-averaged vertical eddy viscosity from the forcing layer.",
            "coarsening_diffusivity_m2_s is the lateral diffusivity used by the merger clock.",
        ],
    }
    for model in ("cl", "scaling"):
        rows = candidate_rows.loc[candidate_rows["model"] == model].copy()
        payload["models"][model] = {
            "n_cases": int(len(rows)),
            "forcing_nu_T_m2_s": _finite_summary(rows["forcing_nu_T_m2_s"]),
            "coarsening_diffusivity_m2_s": _finite_summary(
                rows["coarsening_diffusivity_m2_s"]
            ),
            "coarsening_anisotropy_ratio": _finite_summary(
                rows["coarsening_anisotropy_ratio"]
            ),
            "tau_coarsening_s": _finite_summary(rows["tau_coarsening_s"]),
            "tau_coarsening_h": _finite_summary(rows["tau_coarsening_s"] / 3600.0),
            "pattern_lifetime_h": _finite_summary(rows["pattern_lifetime_s"] / 3600.0),
            "tau_initial_le_lifetime_cases": int(
                (
                    rows["tau_coarsening_s"]
                    <= rows["pattern_lifetime_s"] + 1.0e-9
                ).sum()
            ),
            "n_events_counts": {
                str(int(index)): int(value)
                for index, value in rows["n_coarsening_events"]
                .value_counts()
                .sort_index()
                .items()
            },
            "visible_limited_cases": int(
                (rows["n_visible_events"] < rows["n_coarsening_events"]).sum()
            ),
            "methods": sorted(
                {
                    str(method)
                    for method in rows["coarsening_diffusivity_method"].dropna().unique()
                }
            ),
        }
    return payload


def _write_wp06_output_package(
    output_path: Path,
    *,
    result: dict[str, Any],
    cases: list[ProvisionalObservationCase | MatchedObservationCase],
    candidate_diagnostics_by_case: dict[str, dict[str, Any]],
    representative_case_ids: list[str],
) -> list[str]:
    """Write the full WP-06 comparison output package."""
    written_paths: list[str] = []

    for model in ("cl", "scaling"):
        written_paths.append(
            str(
                write_predicted_vs_observed_plot(
                    result["prediction_table"],
                    model,
                    output_path / f"predicted_vs_observed_{model}.png",
                )
            )
        )
        written_paths.append(
            str(
                write_spacing_vs_wind_plot(
                    result["prediction_table"],
                    model,
                    output_path / f"spacing_vs_wind_{model}.png",
                )
            )
        )
        written_paths.append(
            str(
                write_attractor_diagnostic(
                    result["prediction_table"],
                    model,
                    output_path / f"attractor_diagnostic_{model}.png",
                )
            )
        )

    written_paths.append(
        str(
            write_dynamic_range_comparison(
                result["metrics_table"],
                output_path / "dynamic_range_comparison.png",
            )
        )
    )
    written_paths.append(
        str(
            write_tail_coverage_comparison(
                result["metrics_table"],
                output_path / "tail_coverage_comparison.png",
            )
        )
    )

    with (output_path / "coarsening_closure_audit.json").open(
        "w",
        encoding="utf-8",
    ) as fh:
        json.dump(
            _coarsening_closure_audit_payload(result["prediction_table"]),
            fh,
            indent=2,
            default=_json_ready,
        )
    written_paths.append(str(output_path / "coarsening_closure_audit.json"))

    case_lookup = {case.observation_id: case for case in cases}
    case_table = result["case_table"].set_index("observation_id", drop=False)
    case_diagnostics_dir = output_path / "case_diagnostics"
    case_diagnostics_dir.mkdir(parents=True, exist_ok=True)

    for observation_id in case_table.index:
        model_rows = result["prediction_table"].loc[
            result["prediction_table"]["observation_id"] == observation_id
        ].copy()
        payload = _case_diagnostic_payload(
            case_table.loc[observation_id],
            model_rows,
            candidate_diagnostics_by_case.get(str(observation_id), {}),
            result["baseline_coefficients"],
            result["scenario"],
            forcing_mode=result["forcing_mode"],
            provisional=bool(result["production_candidate_assessment"]["provisional"]),
            representative_case_ids=representative_case_ids,
        )
        case_file = case_diagnostics_dir / f"case_{observation_id}.json"
        with case_file.open("w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, default=_json_ready)
        written_paths.append(str(case_file))

    for observation_id in representative_case_ids:
        timeline = _representative_case_timeline(
            case_lookup[observation_id],
            onset_only=bool(result["scenario"]["onset_only"]),
            visible_spacing_multiplier=float(result["scenario"]["visible_spacing_multiplier"]),
            max_visible_mergers=int(result["scenario"]["max_visible_mergers"]),
            max_visible_aspect_ratio=float(result["scenario"]["max_visible_aspect_ratio"]),
            max_cell_aspect_ratio=float(result["scenario"]["max_cell_aspect_ratio"]),
            drag_method=str(result["scenario"]["drag_method"]),
            drift_method=str(result["scenario"]["drift_method"]),
        )
        written_paths.append(
            str(
                write_enhancement_index_timeseries(
                    timeline,
                    observation_id,
                    output_path / f"enhancement_index_timeseries_{observation_id}.png",
                )
            )
        )

    with (output_path / "production_candidate_assessment.json").open(
        "w",
        encoding="utf-8",
    ) as fh:
        json.dump(
            {
                "production_candidate_assessment": result["production_candidate_assessment"],
                "production_reproducibility": result["production_reproducibility"],
                "assumptions": _comparison_assumptions_payload(),
            },
            fh,
            indent=2,
            default=_json_ready,
        )
    written_paths.append(str(output_path / "production_candidate_assessment.json"))
    return written_paths


def run_wp05_comparison(
    observations: str | Path | pd.DataFrame = "data/raw/observations.csv",
    output_dir: str | Path | None = None,
    site_registry: dict[str, SiteProxySpecification] | None = None,
    pattern_lifetime_s: float = 1800.0,
    era5_cache_dir: str | Path | None = None,
    lookback_hours: float = 6.0,
    onset_only: bool = False,
    visible_spacing_multiplier: float = 1.0,
    max_visible_mergers: int = 3,
    max_visible_aspect_ratio: float = 300.0,
    max_cell_aspect_ratio: float = 12.0,
    drag_method: str = "coare35",
    drift_method: str = "webb_fox_kemper",
) -> dict[str, Any]:
    """
    Run the remaining WP-05 full observation-set comparison.

    Parameters:
        observations:        Benchmark observation CSV or dataframe
        output_dir:          Optional directory for CSV/JSON outputs
        site_registry:       Optional proxy site registry override
        pattern_lifetime_s:  Fixed organisation time for the proxy runner [s]
        era5_cache_dir:      Optional cache directory for forcing-matched reruns
        lookback_hours:      History window used by the matched rerun [h]
        onset_only:          If True, the CL candidate reports raw onset width only
        visible_spacing_multiplier: Common observation-scale spacing multiplier [-]
        max_visible_mergers:        Maximum visible merger levels carried [-]
        max_visible_aspect_ratio:   Observation-scale safety cap [width/depth]
        max_cell_aspect_ratio:      Mechanical cell-scale cap [width/depth]
        drag_method:                Forcing drag parameterisation
        drift_method:               Forcing Stokes-drift parameterisation

    Returns:
        Dict containing case inputs, prediction table, metrics table, baseline
        coefficients, candidate summaries, and the decision gate.
    """
    observation_frame = load_observations(observations)
    if era5_cache_dir is None:
        forcing_mode = "provisional_site_proxy"
        provisional = True
        cases = [
            build_provisional_case(
                observation_row=row,
                site_registry=site_registry,
                pattern_lifetime_s=pattern_lifetime_s,
            )
            for _, row in observation_frame.iterrows()
        ]
    else:
        forcing_mode = "era5_cache_matched"
        provisional = False
        cases = [
            build_matched_case(
                observation_row=row,
                cache_dir=era5_cache_dir,
                site_registry=site_registry,
                lookback_hours=lookback_hours,
            )
            for _, row in observation_frame.iterrows()
        ]
    case_records = []
    candidate_rows = []
    candidate_diagnostics_by_case: dict[str, dict[str, Any]] = {}
    for case, (_, row) in zip(cases, observation_frame.iterrows()):
        case_records.append(
            {
                "observation_id": case.observation_id,
                "site_id": case.site_id,
                "site_label": case.site_label,
                "image_date": case.image_date,
                "observed_spacing_m": case.observed_spacing_m,
                "measurement_method": case.measurement_method,
                "manual_spacing_m": row["manual_spacing_m"],
                "wiggle_spacing_m": row["wiggle_spacing_m"],
                "has_manual": bool(row["has_manual"]),
                "has_wiggle": bool(row["has_wiggle"]),
                "latitude": case.latitude,
                "longitude": case.longitude,
                "observation_time_utc": case.observation_time_utc,
                "representative_u10_mps": case.representative_u10_mps,
                "depth_m": case.depth_m,
                "fetch_m": case.fetch_m,
                "pattern_lifetime_s": case.pattern_lifetime_s,
                "pattern_reset_time": float("nan"),
                "pattern_reset_event_count": 0,
                "pattern_reset_causes": "",
                "wind_direction_deg": case.wind_direction_deg,
                "forcing_mode": case.forcing_mode,
                "matched_cache_file": getattr(case, "cache_file", ""),
                "lookback_hours": getattr(case, "lookback_hours", float("nan")),
                "onset_only": bool(onset_only),
                "visible_spacing_multiplier": float(visible_spacing_multiplier),
                "max_visible_mergers": int(max_visible_mergers),
                "max_visible_aspect_ratio": float(max_visible_aspect_ratio),
                "max_cell_aspect_ratio": float(max_cell_aspect_ratio),
                "drag_method": drag_method,
                "drift_method": drift_method,
                "proxy_notes": case.proxy_notes,
            }
        )
        case_diagnostics: dict[str, Any] = {}
        cl_result = _run_candidate_model(
            case,
            "cl",
            onset_only=onset_only,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
            drag_method=drag_method,
            drift_method=drift_method,
        )
        case_diagnostics["cl"] = cl_result
        case_records[-1]["pattern_lifetime_s"] = float(
            cl_result.get("pattern_lifetime_s", case.pattern_lifetime_s)
        )
        disruption = (
            cl_result.get("disruption", {})
            if isinstance(cl_result.get("disruption"), dict)
            else {}
        )
        reset_causes = disruption.get("reset_causes", [])
        case_records[-1]["pattern_reset_time"] = disruption.get("reset_time", float("nan"))
        case_records[-1]["pattern_reset_event_count"] = int(
            disruption.get("event_count", 0)
        )
        case_records[-1]["pattern_reset_causes"] = (
            "|".join(str(item) for item in reset_causes)
            if isinstance(reset_causes, list)
            else str(reset_causes)
        )
        candidate_rows.append(
            _candidate_prediction_row(
                case,
                row,
                model="cl",
                onset_only=onset_only,
                visible_spacing_multiplier=visible_spacing_multiplier,
                max_visible_mergers=max_visible_mergers,
                max_visible_aspect_ratio=max_visible_aspect_ratio,
                max_cell_aspect_ratio=max_cell_aspect_ratio,
                drag_method=drag_method,
                drift_method=drift_method,
                result=cl_result,
            )
        )
        scaling_result = _run_candidate_model(
            case,
            "scaling",
            onset_only=False,
            visible_spacing_multiplier=visible_spacing_multiplier,
            max_visible_mergers=max_visible_mergers,
            max_visible_aspect_ratio=max_visible_aspect_ratio,
            max_cell_aspect_ratio=max_cell_aspect_ratio,
            drag_method=drag_method,
            drift_method=drift_method,
        )
        case_diagnostics["scaling"] = scaling_result
        candidate_rows.append(
            _candidate_prediction_row(
                case,
                row,
                model="scaling",
                onset_only=False,
                visible_spacing_multiplier=visible_spacing_multiplier,
                max_visible_mergers=max_visible_mergers,
                max_visible_aspect_ratio=max_visible_aspect_ratio,
                max_cell_aspect_ratio=max_cell_aspect_ratio,
                drag_method=drag_method,
                drift_method=drift_method,
                result=scaling_result,
            )
        )
        candidate_diagnostics_by_case[str(case.observation_id)] = case_diagnostics

    case_table = pd.DataFrame(case_records).set_index("observation_id", drop=False)
    candidate_table = pd.DataFrame(candidate_rows)
    baseline_rows, baseline_fits = _baseline_prediction_rows(case_table)
    prediction_table = pd.concat(
        [candidate_table, pd.DataFrame(baseline_rows)],
        ignore_index=True,
    )

    metrics_rows = []
    subsets = build_benchmark_subsets(case_table)
    for subset_id, observation_ids in subsets.items():
        if len(observation_ids) == 0:
            continue
        observed_subset = case_table.loc[observation_ids, "observed_spacing_m"].to_numpy(
            dtype=float
        )
        for model in prediction_table["model"].unique():
            model_rows = (
                prediction_table.loc[prediction_table["model"] == model]
                .set_index("observation_id", drop=False)
                .loc[observation_ids]
            )
            metrics_rows.append(
                _metrics_row(
                    model=model,
                    model_type=str(model_rows["model_type"].iloc[0]),
                    subset_id=subset_id,
                    predicted=model_rows["comparison_spacing_m"].to_numpy(dtype=float),
                    observed=observed_subset,
                )
            )
    metrics_table = pd.DataFrame(metrics_rows).sort_values(
        ["subset_id", "model_type", "model"]
    )

    candidate_summaries = {
        model: _candidate_summary(metrics_table, prediction_table, model)
        for model in ("cl", "scaling")
    }
    decision_gate = _decision_gate_summary(
        metrics_table,
        prediction_table,
        forcing_mode=forcing_mode,
        provisional=provisional,
    )
    production_candidate_assessment = _assess_production_candidate(
        metrics_table,
        forcing_mode=forcing_mode,
        provisional=provisional,
    )
    representative_case_ids = _select_representative_case_ids(case_table)

    result = {
        "forcing_mode": forcing_mode,
        "scenario": {
            "onset_only": bool(onset_only),
            "visible_spacing_multiplier": float(visible_spacing_multiplier),
            "max_visible_mergers": int(max_visible_mergers),
            "max_visible_aspect_ratio": float(max_visible_aspect_ratio),
            "max_cell_aspect_ratio": float(max_cell_aspect_ratio),
            "drag_method": drag_method,
            "drift_method": drift_method,
            "lookback_hours": float(lookback_hours),
        },
        "cases": cases,
        "case_table": case_table.reset_index(drop=True),
        "prediction_table": prediction_table,
        "metrics_table": metrics_table.reset_index(drop=True),
        "baseline_coefficients": {
            model: fit["coefficients"] for model, fit in baseline_fits.items()
        },
        "candidate_summaries": candidate_summaries,
        "decision_gate": decision_gate,
        "production_candidate_assessment": production_candidate_assessment,
        "representative_case_ids": representative_case_ids,
    }
    if era5_cache_dir is not None:
        result["era5_cache_dir"] = str(era5_cache_dir)

    selected_candidate = production_candidate_assessment["selected_candidate"]
    if selected_candidate is None:
        result["production_reproducibility"] = {
            "status": "skipped",
            "candidate": None,
            "reason": production_candidate_assessment["promotion_blocker"],
        }
    else:
        result["production_reproducibility"] = verify_production_candidate_reproduction(
            result,
            candidate=str(selected_candidate),
        )

    if output_dir is not None:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        result["case_table"].to_csv(output_path / "case_inputs.csv", index=False)
        result["prediction_table"].to_csv(
            output_path / "observation_predictions.csv",
            index=False,
        )
        result["metrics_table"].to_csv(output_path / "metrics_table.csv", index=False)
        with (output_path / "decision_gate_summary.json").open("w", encoding="utf-8") as fh:
            json.dump(
                {
                    "decision_gate": decision_gate,
                    "candidate_summaries": candidate_summaries,
                    "baseline_coefficients": result["baseline_coefficients"],
                    "scenario": result["scenario"],
                },
                fh,
                indent=2,
                default=_json_ready,
            )
        result["written_output_paths"] = _write_wp06_output_package(
            output_path,
            result=result,
            cases=cases,
            candidate_diagnostics_by_case=candidate_diagnostics_by_case,
            representative_case_ids=representative_case_ids,
        )
        result["output_dir"] = str(output_path)

    return result
