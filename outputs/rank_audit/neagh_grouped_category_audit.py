#!/usr/bin/env python3
"""
Grouped-point Neagh category audit using point_summary_enriched.csv.

Assumptions made explicit in this script:
- `bathymetry` in point_summary_enriched.csv is stored as negative depth below
  the water surface, so the physical model depth is `abs(bathymetry)` [m].
- `fetch_0` ... `fetch_340` are radial fetch distances indexed by compass
  bearing in 20 degree increments, defined relative to the matched point.
  Example: `fetch_0` is the over-water path from the point toward geographic
  north, so it is the relevant fetch bin for a northerly wind.
- ERA5/Open-Meteo `wind_direction_10m` is treated as the meteorological
  direction the wind is coming from, so the fetch lookup uses the same
  bearing by default.
- `dist_from_prev_m` is treated as the grouped observation-scale spacing for a
  `group_id`, i.e. the grouped mean spacing retained in point_summary_enriched.
- Rows without `observation_date` are excluded because the weather cache key
  cannot be resolved for them.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import neagh_directional_rerun as ndr


CATEGORY_ORDER = ["manual", "stream", "wiggle"]
CATEGORY_COLORS = {
    "manual": "#1f5aa6",
    "stream": "#c97900",
    "wiggle": "#0a7b55",
}
CATEGORY_MARKERS = {
    "manual": "o",
    "stream": "s",
    "wiggle": "^",
}
LAYER_COLUMNS = ["L_inst", "coarsened_width", "mechanical_cell_width", "comparison_spacing"]
LAYER_LABELS = {
    "L_inst": "Onset width L_inst",
    "coarsened_width": "Raw coarsened width",
    "mechanical_cell_width": "Mechanical capped width",
    "comparison_spacing": "Comparison spacing",
}


def _time_weighted_hours(
    timestamps: pd.Series,
    mask: pd.Series,
) -> float:
    """Return the time spent in rows selected by `mask` [h]."""
    times = pd.to_datetime(timestamps, utc=True)
    deltas = (times.shift(-1) - times).dt.total_seconds().fillna(0.0) / 3600.0
    return float(deltas.loc[mask].sum())


def load_grouped_points(point_summary_path: Path) -> pd.DataFrame:
    """Load grouped Neagh points with valid dates and observed spacing."""
    frame = pd.read_csv(point_summary_path)
    frame = frame.loc[frame["observation_date"].notna()].copy()
    frame["observation_date"] = pd.to_datetime(frame["observation_date"])
    frame["observed_spacing"] = frame["dist_from_prev_m"].astype(float)
    frame["depth"] = frame["bathymetry"].astype(float).abs()
    frame["site_id"] = "neagh"
    return frame.sort_values(["observation_date", "category", "group_id"]).reset_index(
        drop=True
    )


def build_point_observation_row(point_row: pd.Series) -> pd.Series:
    """Convert one grouped point row to the observation schema used by the rerun."""
    category = str(point_row["category"])
    observed_spacing = float(point_row["observed_spacing"])
    data: dict[str, Any] = {
        "observation_id": str(point_row["group_id"]),
        "image_date": pd.Timestamp(point_row["observation_date"]),
        "authoritative_lat": float(point_row["lat"]),
        "authoritative_lng": float(point_row["lng"]),
        "target_spacing_m": observed_spacing,
        "measurement_method": category,
        "manual_spacing_m": observed_spacing if category in {"manual", "stream"} else float("nan"),
        "wiggle_spacing_m": observed_spacing if category == "wiggle" else float("nan"),
        "point_group_id": str(point_row["group_id"]),
        "annotation_category": category,
        "point_match_distance_m": 0.0,
        "point_length_m": ndr._as_float(point_row.get("length_m")),
        "point_bearing_deg": ndr._as_float(point_row.get("bearing_deg")),
        "point_mean_spacing_m": observed_spacing,
        "point_group_index_mean": ndr._as_float(point_row.get("group_index")),
        "point_n_annotations": ndr._as_float(point_row.get("n_points")),
        "bathymetry_raw_m": float(point_row["bathymetry"]),
        "depth_m": float(point_row["depth"]),
        "observation_label": str(point_row.get("label", "")),
        "first_timestamp": str(point_row.get("first_timestamp", "")),
        "last_timestamp": str(point_row.get("last_timestamp", "")),
        "source_files": str(point_row.get("source_files", "")),
    }
    for bearing in range(0, 360, 20):
        data[f"fetch_{bearing}"] = float(point_row[f"fetch_{bearing}"])
    return pd.Series(data)


def rerun_grouped_point(
    point_row: pd.Series,
    cache_dir: Path,
    lookback_hours: float,
    angle_mode: str,
    drag_method: str,
    drift_method: str,
) -> dict[str, Any]:
    """Run the CL analysis for one grouped point row with directional fetch."""
    observation_row = build_point_observation_row(point_row)
    registry = ndr.default_site_proxy_registry()
    neagh_spec = registry["neagh"]
    observation_time_utc = ndr.observation_time_from_proxy(
        image_date=pd.Timestamp(observation_row["image_date"]),
        site_spec=neagh_spec,
    )
    cache_file = ndr._cache_file_for_observation(observation_row, cache_dir)
    if not cache_file.exists():
        raise FileNotFoundError(
            f"ERA5 cache file missing for grouped point {observation_row['observation_id']}: "
            f"{cache_file}"
        )

    era5_frame = ndr._load_era5_cache_frame(cache_file)
    weather_data = ndr._weather_history_with_observation_sample(
        era5_frame=era5_frame,
        observation_time_utc=observation_time_utc,
        lookback_hours=lookback_hours,
    )
    forcing_history, disruption_history, history_frame = ndr._forcing_history_with_directional_fetch(
        weather_data=weather_data,
        observation_time_utc=observation_time_utc,
        depth_m=float(observation_row["depth_m"]),
        drag_method=drag_method,
        drift_method=drift_method,
        angle_mode=angle_mode,
        observation_row=observation_row,
        lookback_hours=lookback_hours,
    )
    pattern_lifetime = ndr._estimate_pattern_lifetime(
        forcing_history=forcing_history,
        disruption_history=disruption_history,
        lookback_hours=lookback_hours,
    )
    disruption = ndr.disruption_check(
        forcing_history=disruption_history,
        lookback_hours=lookback_hours,
    )
    result = ndr.analyse_candidate_cl(
        forcing=forcing_history[-1],
        pattern_lifetime=pattern_lifetime,
        environmental=ndr.build_environmental_context(),
        onset_only=False,
        visible_spacing_multiplier=1.0,
        max_visible_mergers=3,
        max_visible_aspect_ratio=300.0,
        max_cell_aspect_ratio=12.0,
    )

    coarsening = result.get("coarsening", {})
    critical_result = result.get("intermediate", {}).get("critical_result", {})
    l_c = ndr._first_finite(
        ndr._field(critical_result, "l_c"),
        ndr._field(critical_result, "lcNL"),
    )
    depth_m = float(observation_row["depth_m"])
    l_inst = (
        float((2.0 * math.pi * depth_m) / l_c)
        if math.isfinite(l_c) and l_c > 0.0
        else ndr._first_finite(coarsening.get("initial_cell_width_m"))
    )
    current_row = history_frame.iloc[-1]
    history_times = pd.to_datetime(history_frame["timestamp"], utc=True)
    history_hours = (
        float((history_times.iloc[-1] - history_times.iloc[0]).total_seconds() / 3600.0)
        if len(history_times) >= 2
        else 0.0
    )
    low4_mask = history_frame["U10"].astype(float) <= 4.0
    low5_mask = history_frame["U10"].astype(float) <= 5.0
    high6_mask = history_frame["U10"].astype(float) >= 6.0
    raw_coarsened_width = ndr._first_finite(coarsening.get("raw_coarsened_width_m"))
    mechanical_cell_width = ndr._first_finite(coarsening.get("coarsened_width_m"))
    comparison_spacing = ndr._first_finite(coarsening.get("visible_spacing_m"))
    visible_hierarchy_width = ndr._first_finite(coarsening.get("visible_hierarchy_width_m"))
    visible_spacing_lower = ndr._first_finite(coarsening.get("visible_spacing_lower_m"))

    return {
        "case_id": str(observation_row["observation_id"]),
        "image_date": str(pd.Timestamp(observation_row["image_date"]).date()),
        "site_id": "neagh",
        "annotation_category": str(observation_row["annotation_category"]),
        "measurement_method": str(observation_row["measurement_method"]),
        "observed_spacing": float(observation_row["target_spacing_m"]),
        "point_group_id": str(observation_row["point_group_id"]),
        "point_n_annotations": ndr._as_float(observation_row["point_n_annotations"]),
        "point_length_m": ndr._as_float(observation_row["point_length_m"]),
        "point_bearing_deg": ndr._as_float(observation_row["point_bearing_deg"]),
        "point_group_index_mean": ndr._as_float(observation_row["point_group_index_mean"]),
        "bathymetry_raw_m": float(observation_row["bathymetry_raw_m"]),
        "depth": depth_m,
        "authoritative_lat": float(observation_row["authoritative_lat"]),
        "authoritative_lng": float(observation_row["authoritative_lng"]),
        "U10": float(result["forcing_summary"]["U10"]),
        "window_mean_u10": float(history_frame["U10"].mean()),
        "window_median_u10": float(history_frame["U10"].median()),
        "window_min_u10": float(history_frame["U10"].min()),
        "window_max_u10": float(history_frame["U10"].max()),
        "window_std_u10": float(history_frame["U10"].std(ddof=0)),
        "window_history_h": history_hours,
        "hours_u10_le_4": _time_weighted_hours(history_frame["timestamp"], low4_mask),
        "hours_u10_le_5": _time_weighted_hours(history_frame["timestamp"], low5_mask),
        "hours_u10_ge_6": _time_weighted_hours(history_frame["timestamp"], high6_mask),
        "wind_direction_deg": float(current_row["wind_direction_deg"]),
        "fetch_bin_deg": int(current_row["fetch_bin_deg"]),
        "fetch_at_observation_m": float(current_row["fetch_m"]),
        "fetch_mean_window_m": float(history_frame["fetch_m"].mean()),
        "fetch_min_window_m": float(history_frame["fetch_m"].min()),
        "fetch_max_window_m": float(history_frame["fetch_m"].max()),
        "Ra": float(result["Ra"]),
        "L_inst": float(l_inst),
        "coarsened_width": raw_coarsened_width,
        "mechanical_cell_width": mechanical_cell_width,
        "visible_hierarchy_width": visible_hierarchy_width,
        "visible_spacing_lower": visible_spacing_lower,
        "comparison_spacing": comparison_spacing,
        "pattern_lifetime_s": float(pattern_lifetime),
        "pattern_lifetime_h": float(pattern_lifetime / 3600.0),
        "u_star_water": float(result["forcing_summary"]["u_star_water"]),
        "nu_T_vertical": float(result["forcing_summary"]["nu_T"]),
        "A_H": float(coarsening.get("coarsening_diffusivity_m2_s", float("nan"))),
        "tau_coarsening_s": float(coarsening.get("tau_coarsening_s", float("nan"))),
        "tau_next_coarsening_s": float(coarsening.get("tau_next_coarsening_s", float("nan"))),
        "n_events": int(coarsening.get("n_events", 0)),
        "visible_n_events": int(coarsening.get("visible_n_events", 0)),
        "cell_cap_binding": bool(coarsening.get("cap_binding", False)),
        "visible_cap_binding": bool(coarsening.get("visible_cap_binding", False)),
        "coarsening_gain_m": float(raw_coarsened_width - l_inst),
        "mechanical_cap_loss_m": float(raw_coarsened_width - mechanical_cell_width),
        "comparison_minus_raw_coarsened_m": float(comparison_spacing - raw_coarsened_width),
        "comparison_minus_mechanical_m": float(comparison_spacing - mechanical_cell_width),
        "comparison_minus_onset_m": float(comparison_spacing - l_inst),
        "regime": str(result["regime"]),
        "disruption_triggered": bool(disruption.get("is_disrupted", False)),
        "cache_file": str(cache_file),
        "source_files": str(observation_row.get("source_files", "")),
    }


def build_layer_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarise observed-vs-layer agreement by annotation category."""
    rows: list[dict[str, Any]] = []
    for category in ["all"] + CATEGORY_ORDER:
        subset = frame if category == "all" else frame.loc[frame["annotation_category"] == category]
        for layer in LAYER_COLUMNS:
            valid = subset[["observed_spacing", layer]].dropna()
            if len(valid) < 2 or valid["observed_spacing"].nunique() < 2 or valid[layer].nunique() < 2:
                rho = float("nan")
                p_value = float("nan")
            else:
                rho, p_value = spearmanr(valid["observed_spacing"], valid[layer])
            bias = valid[layer] - valid["observed_spacing"]
            ratio = valid[layer] / valid["observed_spacing"]
            rows.append(
                {
                    "annotation_category": category,
                    "layer": layer,
                    "n": int(len(valid)),
                    "spearman_rho": ndr._as_float(rho),
                    "spearman_p_value": ndr._as_float(p_value),
                    "median_observed_m": float(valid["observed_spacing"].median()) if not valid.empty else float("nan"),
                    "median_layer_m": float(valid[layer].median()) if not valid.empty else float("nan"),
                    "median_bias_m": float(bias.median()) if not valid.empty else float("nan"),
                    "mae_m": float(bias.abs().mean()) if not valid.empty else float("nan"),
                    "median_ratio": float(ratio.median()) if not valid.empty else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def build_wind_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarise wind-duration/coarsening diagnostics by annotation category."""
    rows: list[dict[str, Any]] = []
    for category in ["all"] + CATEGORY_ORDER:
        subset = frame if category == "all" else frame.loc[frame["annotation_category"] == category]
        if subset.empty:
            continue
        rows.append(
            {
                "annotation_category": category,
                "n": int(len(subset)),
                "median_observed_m": float(subset["observed_spacing"].median()),
                "median_U10_mps": float(subset["U10"].median()),
                "median_window_mean_u10_mps": float(subset["window_mean_u10"].median()),
                "median_pattern_lifetime_h": float(subset["pattern_lifetime_h"].median()),
                "median_hours_u10_le_4": float(subset["hours_u10_le_4"].median()),
                "median_hours_u10_le_5": float(subset["hours_u10_le_5"].median()),
                "median_hours_u10_ge_6": float(subset["hours_u10_ge_6"].median()),
                "median_L_inst_m": float(subset["L_inst"].median()),
                "median_coarsened_width_m": float(subset["coarsened_width"].median()),
                "median_comparison_spacing_m": float(subset["comparison_spacing"].median()),
                "median_n_events": float(subset["n_events"].median()),
                "fraction_n_events_zero": float((subset["n_events"] == 0).mean()),
                "median_coarsening_gain_m": float(subset["coarsening_gain_m"].median()),
                "median_comparison_minus_raw_coarsened_m": float(
                    subset["comparison_minus_raw_coarsened_m"].median()
                ),
                "median_comparison_minus_mechanical_m": float(
                    subset["comparison_minus_mechanical_m"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def build_best_layer_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """Count which layer is closest to the grouped observed spacing."""
    enriched = frame.copy()
    error_cols = []
    for layer in LAYER_COLUMNS:
        col = f"abs_error_{layer}"
        enriched[col] = (enriched[layer] - enriched["observed_spacing"]).abs()
        error_cols.append(col)
    best = enriched[error_cols].idxmin(axis=1).str.replace("abs_error_", "", regex=False)
    enriched["closest_layer"] = best
    rows: list[dict[str, Any]] = []
    for category in ["all"] + CATEGORY_ORDER:
        subset = enriched if category == "all" else enriched.loc[enriched["annotation_category"] == category]
        if subset.empty:
            continue
        counts = subset["closest_layer"].value_counts()
        for layer in LAYER_COLUMNS:
            rows.append(
                {
                    "annotation_category": category,
                    "closest_layer": layer,
                    "count": int(counts.get(layer, 0)),
                }
            )
    return pd.DataFrame(rows)


def _annotate_small_group(ax: plt.Axes, subset: pd.DataFrame) -> None:
    """Label small groups directly to keep wiggle cases visible."""
    if len(subset) > 5:
        return
    for _, row in subset.iterrows():
        ax.annotate(
            str(row["image_date"]),
            (float(row["_x"]), float(row["_y"])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
            color=CATEGORY_COLORS[str(row["annotation_category"])],
        )


def plot_observed_vs_layers(frame: pd.DataFrame, output_path: Path) -> None:
    """Scatter observed spacing against each model layer by annotation category."""
    fig, axes = plt.subplots(1, len(LAYER_COLUMNS), figsize=(20.0, 4.8))
    for ax, layer in zip(axes, LAYER_COLUMNS, strict=True):
        combined_min = min(frame["observed_spacing"].min(), frame[layer].min())
        combined_max = max(frame["observed_spacing"].max(), frame[layer].max())
        pad = 0.05 * (combined_max - combined_min)
        line = np.array([combined_min - pad, combined_max + pad])
        ax.plot(line, line, linestyle="--", linewidth=1.0, color="#555555")
        for category in CATEGORY_ORDER:
            subset = frame.loc[frame["annotation_category"] == category].copy()
            if subset.empty:
                continue
            subset["_x"] = subset["observed_spacing"]
            subset["_y"] = subset[layer]
            ax.scatter(
                subset["_x"],
                subset["_y"],
                s=46,
                alpha=0.9,
                linewidth=0.5,
                edgecolors="white",
                color=CATEGORY_COLORS[category],
                marker=CATEGORY_MARKERS[category],
                label=category,
            )
            _annotate_small_group(ax, subset)
        ax.set_title(LAYER_LABELS[layer])
        ax.set_xlabel("Observed grouped spacing [m]")
        ax.set_ylabel("Model spacing [m]")
        ax.grid(alpha=0.2, linewidth=0.6)
    handles = [
        Line2D(
            [0],
            [0],
            marker=CATEGORY_MARKERS[category],
            color="none",
            markerfacecolor=CATEGORY_COLORS[category],
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=8,
            label=category,
        )
        for category in CATEGORY_ORDER
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Grouped Neagh observations vs model layers", y=1.02)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_layer_bias_by_category(frame: pd.DataFrame, output_path: Path) -> None:
    """Strip plot of layer bias relative to grouped observed spacing."""
    rng = np.random.default_rng(42)
    fig, axes = plt.subplots(1, len(LAYER_COLUMNS), figsize=(20.0, 4.8), sharey=True)
    positions = {category: idx for idx, category in enumerate(CATEGORY_ORDER)}
    for ax, layer in zip(axes, LAYER_COLUMNS, strict=True):
        for category in CATEGORY_ORDER:
            subset = frame.loc[frame["annotation_category"] == category].copy()
            if subset.empty:
                continue
            bias = subset[layer] - subset["observed_spacing"]
            jitter = rng.uniform(-0.12, 0.12, size=len(subset))
            x = positions[category] + jitter
            ax.scatter(
                x,
                bias,
                s=42,
                alpha=0.85,
                linewidth=0.5,
                edgecolors="white",
                color=CATEGORY_COLORS[category],
                marker=CATEGORY_MARKERS[category],
            )
            median_bias = float(bias.median())
            ax.hlines(
                median_bias,
                positions[category] - 0.22,
                positions[category] + 0.22,
                color="#222222",
                linewidth=2.0,
            )
        ax.axhline(0.0, color="#555555", linewidth=1.0, linestyle="--")
        ax.set_title(LAYER_LABELS[layer])
        ax.set_xticks(list(positions.values()), CATEGORY_ORDER)
        ax.set_xlabel("Annotation category")
        ax.grid(axis="y", alpha=0.2, linewidth=0.6)
    axes[0].set_ylabel("Layer bias = model - observed [m]")
    fig.suptitle("Where the overprediction enters by category", y=1.02)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_wind_duration_diagnostics(frame: pd.DataFrame, output_path: Path) -> None:
    """Plot wind-speed/duration diagnostics against observed spacing and error ratio."""
    plot_frame = frame.copy()
    plot_frame["comparison_ratio"] = plot_frame["comparison_spacing"] / plot_frame["observed_spacing"]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))
    specs = [
        ("window_mean_u10", "observed_spacing", "Window-mean U10 [m/s]", "Observed grouped spacing [m]"),
        ("pattern_lifetime_h", "observed_spacing", "Pattern lifetime [h]", "Observed grouped spacing [m]"),
        ("window_mean_u10", "comparison_ratio", "Window-mean U10 [m/s]", "Comparison / observed [-]"),
        ("pattern_lifetime_h", "comparison_ratio", "Pattern lifetime [h]", "Comparison / observed [-]"),
    ]
    for ax, (x_col, y_col, x_label, y_label) in zip(axes.flat, specs, strict=True):
        for category in CATEGORY_ORDER:
            subset = plot_frame.loc[plot_frame["annotation_category"] == category].copy()
            if subset.empty:
                continue
            subset["_x"] = subset[x_col]
            subset["_y"] = subset[y_col]
            ax.scatter(
                subset["_x"],
                subset["_y"],
                s=54,
                alpha=0.9,
                linewidth=0.5,
                edgecolors="white",
                color=CATEGORY_COLORS[category],
                marker=CATEGORY_MARKERS[category],
            )
            _annotate_small_group(ax, subset)
        if y_col == "comparison_ratio":
            ax.axhline(1.0, color="#555555", linewidth=1.0, linestyle="--")
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid(alpha=0.2, linewidth=0.6)
    handles = [
        Line2D(
            [0],
            [0],
            marker=CATEGORY_MARKERS[category],
            color="none",
            markerfacecolor=CATEGORY_COLORS[category],
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=8,
            label=category,
        )
        for category in CATEGORY_ORDER
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Wind-speed and duration diagnostics by category", y=0.995)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def print_summary(
    frame: pd.DataFrame,
    layer_summary: pd.DataFrame,
    wind_summary: pd.DataFrame,
    best_layer_summary: pd.DataFrame,
    excluded_missing_date: int,
) -> None:
    """Print the high-signal grouped-point audit summary."""
    print(f"Grouped Neagh rows rerun: {len(frame)}")
    print(f"Excluded rows with missing observation_date: {excluded_missing_date}")
    print(
        "Category counts: "
        + ", ".join(
            f"{category}={int((frame['annotation_category'] == category).sum())}"
            for category in CATEGORY_ORDER
        )
    )
    print(
        "Observation-scale spacing [m]: "
        f"min={frame['observed_spacing'].min():.1f}, "
        f"median={frame['observed_spacing'].median():.1f}, "
        f"max={frame['observed_spacing'].max():.1f}"
    )
    print(
        "Comparison spacing - raw coarsened width [m]: "
        f"median={frame['comparison_minus_raw_coarsened_m'].median():.1f}, "
        f"max={frame['comparison_minus_raw_coarsened_m'].max():.1f}"
    )
    print(
        "Comparison spacing - mechanical capped width [m]: "
        f"median={frame['comparison_minus_mechanical_m'].median():.1f}, "
        f"max={frame['comparison_minus_mechanical_m'].max():.1f}"
    )
    print("\nLayer summary")
    print(layer_summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print("\nWind/duration summary")
    print(wind_summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print("\nClosest layer counts")
    print(best_layer_summary.to_string(index=False))


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--point-summary",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "point_summary_enriched.csv",
        help="Grouped point summary with bathymetry and directional fetch [path]",
    )
    parser.add_argument(
        "--era5-cache-dir",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "era5_cache",
        help="Directory containing grouped-point ERA5/Open-Meteo cache JSONs [path]",
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
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs" / "rank_audit" / "neagh_grouped_category_audit",
        help="Directory for the grouped-point audit tables and plots [path]",
    )
    return parser.parse_args()


def main() -> None:
    """Run the grouped-point category audit and write tables/plots."""
    args = parse_args()
    points = load_grouped_points(args.point_summary)
    total_rows = len(pd.read_csv(args.point_summary))
    excluded_missing_date = int(total_rows - len(points))

    rows = [
        rerun_grouped_point(
            point_row=row,
            cache_dir=args.era5_cache_dir,
            lookback_hours=args.lookback_hours,
            angle_mode=args.fetch_angle_mode,
            drag_method="coare35",
            drift_method="webb_fox_kemper",
        )
        for _, row in points.iterrows()
    ]
    frame = pd.DataFrame(rows).sort_values(["annotation_category", "image_date", "case_id"]).reset_index(drop=True)
    frame["annotation_category"] = pd.Categorical(
        frame["annotation_category"],
        categories=CATEGORY_ORDER,
        ordered=True,
    )

    layer_summary = build_layer_summary(frame)
    wind_summary = build_wind_summary(frame)
    best_layer_summary = build_best_layer_summary(frame)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "neagh_grouped_category_rank_audit.csv", index=False)
    layer_summary.to_csv(args.output_dir / "neagh_grouped_category_layer_summary.csv", index=False)
    wind_summary.to_csv(args.output_dir / "neagh_grouped_category_wind_summary.csv", index=False)
    best_layer_summary.to_csv(
        args.output_dir / "neagh_grouped_category_best_layer_summary.csv",
        index=False,
    )

    plot_observed_vs_layers(
        frame=frame,
        output_path=args.output_dir / "observed_vs_layers_by_category.png",
    )
    plot_layer_bias_by_category(
        frame=frame,
        output_path=args.output_dir / "layer_bias_by_category.png",
    )
    plot_wind_duration_diagnostics(
        frame=frame,
        output_path=args.output_dir / "wind_duration_diagnostics.png",
    )

    print(f"Wrote grouped audit to {args.output_dir}")
    print_summary(
        frame=frame,
        layer_summary=layer_summary,
        wind_summary=wind_summary,
        best_layer_summary=best_layer_summary,
        excluded_missing_date=excluded_missing_date,
    )


if __name__ == "__main__":
    main()
