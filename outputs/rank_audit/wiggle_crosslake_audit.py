#!/usr/bin/env python3
"""
Cross-lake wiggle audit using benchmark observations.csv and current matched rerun outputs.

Purpose:
- Use `wiggle_spacing_m` as the explicit observation target whenever it exists,
  even for mixed manual+wiggle rows.
- Test whether the apparent wiggle agreement seen at Neagh also appears at
  Erie, Prairie, and Taihu.
- Add matched wind-window diagnostics and time-ordered plots so the wiggle
  spacing can be compared against the forcing conditions that preceded each
  observation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.evaluation.comparison import (
    build_matched_case,
    default_site_proxy_registry,
    load_observations,
)


SITE_ORDER = ["erie", "neagh", "prairie", "taihu"]
SITE_LABELS = {
    "erie": "Erie",
    "neagh": "Neagh",
    "prairie": "Prairie",
    "taihu": "Taihu",
}
SITE_COLORS = {
    "erie": "#7a3e9d",
    "neagh": "#0f6d5b",
    "prairie": "#b45f06",
    "taihu": "#1d5ea8",
}
LAYER_COLUMNS = [
    "L_inst",
    "mechanical_capped_width",
    "comparison_spacing",
]


def _as_float(value: Any) -> float:
    """Convert scalar JSON-like values to float, preserving NaN on failure."""
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


def _as_int(value: Any, default: int = 0) -> int:
    """Convert a numeric value to int with a fallback."""
    number = _as_float(value)
    if math.isfinite(number):
        return int(number)
    return int(default)


def _as_bool(value: Any, default: bool = False) -> bool:
    """Convert common scalar representations to bool."""
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        if math.isnan(float(value)):
            return bool(default)
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes"}:
            return True
        if lowered in {"false", "0", "no", ""}:
            return False
    return bool(default)


def _resolve_case_dir(path: Path) -> Path:
    """Resolve either a comparison output directory or a direct case_diagnostics dir."""
    if path.name == "case_diagnostics":
        case_dir = path
    else:
        case_dir = path / "case_diagnostics"
    if not case_dir.exists():
        raise FileNotFoundError(f"Case diagnostics directory not found: {case_dir}")
    return case_dir


def _compute_l_inst(depth_m: Any, critical_result: dict[str, Any], onset_width_m: Any) -> float:
    """Compute the onset width L_inst [m] from the critical wavenumber when possible."""
    depth = _as_float(depth_m)
    l_c = _first_finite(critical_result.get("l_c"), critical_result.get("lcNL"))
    if math.isfinite(depth) and depth > 0.0 and math.isfinite(l_c) and l_c > 0.0:
        return float((2.0 * math.pi * depth) / l_c)
    return _as_float(onset_width_m)


def _weighted_hours(timestamps: pd.Series, mask: pd.Series) -> float:
    """Estimate hours spent in rows selected by `mask`."""
    times = pd.to_datetime(timestamps, utc=True)
    deltas = (times.shift(-1) - times).dt.total_seconds().fillna(0.0) / 3600.0
    return float(deltas.loc[mask].sum())


def load_wiggle_observations(observation_path: Path) -> pd.DataFrame:
    """Load only the benchmark rows that contain a wiggle spacing."""
    frame = load_observations(observation_path)
    wiggle = frame.loc[frame["wiggle_spacing_m"].notna()].copy().reset_index(drop=True)
    registry = default_site_proxy_registry()
    wiggle["site_id"] = wiggle.apply(
        lambda row: next(
            (
                spec.site_id
                for spec in registry.values()
                if spec.contains(
                    latitude=float(row["authoritative_lat"]),
                    longitude=float(row["authoritative_lng"]),
                )
            ),
            "unknown",
        ),
        axis=1,
    )
    return wiggle.sort_values(["site_id", "image_date", "observation_id"]).reset_index(
        drop=True
    )


def extract_case_row(
    case_file: Path,
    observation_row: pd.Series,
) -> dict[str, Any]:
    """Flatten one current matched case JSON into the wiggle-audit fields."""
    with case_file.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)

    model_payload = payload["models"]["cl"]
    prediction_row = model_payload.get("prediction_row", {})
    full = model_payload.get("full_diagnostic", {})
    case_input = payload.get("case_input", {})
    forcing = full.get("forcing_summary", {})
    coarsening = full.get("coarsening", {})
    intermediate = full.get("intermediate", {})
    critical_result = intermediate.get("critical_result", {})

    depth = _first_finite(
        forcing.get("depth"),
        case_input.get("depth_m"),
        prediction_row.get("depth_m"),
    )
    onset_width = _first_finite(
        coarsening.get("initial_cell_width_m"),
        prediction_row.get("initial_cell_width_m"),
    )
    n_events = _as_int(
        _first_finite(coarsening.get("n_events"), prediction_row.get("n_coarsening_events"))
    )
    comparison_spacing = _first_finite(
        coarsening.get("visible_spacing_m"),
        prediction_row.get("comparison_spacing_m"),
        prediction_row.get("raw_predicted_spacing_m"),
    )
    return {
        "case_id": str(observation_row["observation_id"]),
        "image_date": pd.Timestamp(observation_row["image_date"]),
        "site_id": str(prediction_row.get("site_id") or case_input.get("site_id") or observation_row["site_id"]),
        "wiggle_spacing": float(observation_row["wiggle_spacing_m"]),
        "manual_spacing": _as_float(observation_row["manual_spacing_m"]),
        "measurement_method": str(observation_row["measurement_method"]),
        "observed_spacing": float(observation_row["wiggle_spacing_m"]),
        "depth": depth,
        "U10": _first_finite(
            forcing.get("U10"),
            case_input.get("representative_u10_mps"),
            prediction_row.get("representative_u10_mps"),
        ),
        "Ra": _first_finite(full.get("Ra"), prediction_row.get("Ra")),
        "L_inst": _compute_l_inst(depth, critical_result, onset_width),
        "coarsened_width": _first_finite(coarsening.get("raw_coarsened_width_m")),
        "mechanical_capped_width": _first_finite(coarsening.get("coarsened_width_m")),
        "comparison_spacing": comparison_spacing,
        "visible_spacing_lower": _first_finite(coarsening.get("visible_spacing_lower_m")),
        "pattern_lifetime_s": _first_finite(
            full.get("pattern_lifetime_s"),
            case_input.get("pattern_lifetime_s"),
            prediction_row.get("pattern_lifetime_s"),
        ),
        "u_star_water": _first_finite(forcing.get("u_star_water")),
        "nu_T_vertical": _first_finite(
            forcing.get("nu_T"),
            coarsening.get("forcing_nu_T_m2_s"),
            prediction_row.get("forcing_nu_T_m2_s"),
        ),
        "A_H": _first_finite(
            coarsening.get("coarsening_diffusivity_m2_s"),
            prediction_row.get("coarsening_diffusivity_m2_s"),
        ),
        "n_events": n_events,
        "visible_n_events": _as_int(
            _first_finite(coarsening.get("visible_n_events"), prediction_row.get("n_visible_events"))
        ),
        "cap_binding": _as_bool(
            coarsening.get("cap_binding", prediction_row.get("cap_binding", False))
        ),
        "regime": str(full.get("regime") or prediction_row.get("regime") or ""),
        "matched_cache_file": str(
            prediction_row.get("matched_cache_file") or case_input.get("cache_file") or ""
        ),
    }


def weather_metrics_for_observation(
    observation_row: pd.Series,
    cache_dir: Path,
    lookback_hours: float,
) -> dict[str, Any]:
    """Compute matched wind-window diagnostics for one observation."""
    matched = build_matched_case(
        observation_row=observation_row,
        cache_dir=cache_dir,
        site_registry=default_site_proxy_registry(),
        lookback_hours=lookback_hours,
    )
    weather = matched.weather_data.copy().sort_values("timestamp").reset_index(drop=True)
    current = weather.iloc[-1]
    low4 = weather["U10"].astype(float) <= 4.0
    low5 = weather["U10"].astype(float) <= 5.0
    high6 = weather["U10"].astype(float) >= 6.0
    history_hours = (
        float(
            (
                pd.to_datetime(weather["timestamp"], utc=True).iloc[-1]
                - pd.to_datetime(weather["timestamp"], utc=True).iloc[0]
            ).total_seconds()
            / 3600.0
        )
        if len(weather) >= 2
        else 0.0
    )
    return {
        "observation_time_utc": matched.observation_time_utc,
        "window_history_h": history_hours,
        "obs_u10": float(current["U10"]),
        "obs_wind_direction_deg": float(current["wind_direction_deg"]),
        "window_mean_u10": float(weather["U10"].mean()),
        "window_median_u10": float(weather["U10"].median()),
        "window_min_u10": float(weather["U10"].min()),
        "window_max_u10": float(weather["U10"].max()),
        "window_std_u10": float(weather["U10"].std(ddof=0)),
        "hours_u10_le_4": _weighted_hours(weather["timestamp"], low4),
        "hours_u10_le_5": _weighted_hours(weather["timestamp"], low5),
        "hours_u10_ge_6": _weighted_hours(weather["timestamp"], high6),
        "pattern_proxy_cache_file": str(matched.cache_file),
    }


def build_wiggle_audit_table(
    case_dir: Path,
    observation_path: Path,
    cache_dir: Path,
    lookback_hours: float,
) -> pd.DataFrame:
    """Build the wiggle-only cross-lake audit table."""
    observations = load_wiggle_observations(observation_path)
    rows: list[dict[str, Any]] = []
    for _, observation_row in observations.iterrows():
        case_file = case_dir / f"case_{observation_row['observation_id']}.json"
        if not case_file.exists():
            raise FileNotFoundError(f"Missing case JSON for {observation_row['observation_id']}: {case_file}")
        row = extract_case_row(case_file=case_file, observation_row=observation_row)
        row.update(
            weather_metrics_for_observation(
                observation_row=observation_row,
                cache_dir=cache_dir,
                lookback_hours=lookback_hours,
            )
        )
        row["has_manual"] = bool(pd.notna(observation_row["manual_spacing_m"]))
        row["pattern_lifetime_h"] = row["pattern_lifetime_s"] / 3600.0
        row["comparison_bias_m"] = row["comparison_spacing"] - row["wiggle_spacing"]
        row["comparison_ratio"] = row["comparison_spacing"] / row["wiggle_spacing"]
        row["onset_bias_m"] = row["L_inst"] - row["wiggle_spacing"]
        row["mechanical_bias_m"] = row["mechanical_capped_width"] - row["wiggle_spacing"]
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["site_id", "image_date", "case_id"]).reset_index(drop=True)


def _spearman(valid: pd.DataFrame, y_col: str, x_col: str) -> tuple[float, float]:
    """Compute Spearman rho/p-value with NaN guards."""
    subset = valid[[y_col, x_col]].dropna()
    if len(subset) < 2 or subset[y_col].nunique() < 2 or subset[x_col].nunique() < 2:
        return float("nan"), float("nan")
    rho, p_value = spearmanr(subset[y_col], subset[x_col])
    return _as_float(rho), _as_float(p_value)


def build_layer_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarise wiggle-target rank agreement overall and by site."""
    rows: list[dict[str, Any]] = []
    summary_sites = ["all"] + SITE_ORDER
    x_cols = LAYER_COLUMNS + [
        "U10",
        "obs_u10",
        "window_mean_u10",
        "pattern_lifetime_h",
        "n_events",
    ]
    for site_id in summary_sites:
        subset = frame if site_id == "all" else frame.loc[frame["site_id"] == site_id]
        for column in x_cols:
            rho, p_value = _spearman(subset, "wiggle_spacing", column)
            valid = subset[["wiggle_spacing", column]].dropna()
            bias = (
                subset[column] - subset["wiggle_spacing"]
                if column in LAYER_COLUMNS
                else pd.Series(dtype=float)
            )
            rows.append(
                {
                    "site_id": site_id,
                    "x": column,
                    "n": int(len(valid)),
                    "rho": rho,
                    "p_value": p_value,
                    "median_wiggle_m": float(subset["wiggle_spacing"].median()) if not subset.empty else float("nan"),
                    "median_x": float(subset[column].median()) if column in subset and not subset.empty else float("nan"),
                    "median_bias_m": float(bias.median()) if not bias.empty else float("nan"),
                    "mae_m": float(bias.abs().mean()) if not bias.empty else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def build_site_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarise wiggle counts, medians, and predicted-vs-observed fit by site."""
    rows: list[dict[str, Any]] = []
    for site_id in SITE_ORDER:
        subset = frame.loc[frame["site_id"] == site_id].copy()
        if subset.empty:
            continue
        rho, p_value = _spearman(subset, "wiggle_spacing", "comparison_spacing")
        rows.append(
            {
                "site_id": site_id,
                "n": int(len(subset)),
                "n_mixed_manual_wiggle": int(subset["has_manual"].sum()),
                "spearman_rho_comparison": rho,
                "spearman_p_value_comparison": p_value,
                "median_wiggle_m": float(subset["wiggle_spacing"].median()),
                "median_comparison_m": float(subset["comparison_spacing"].median()),
                "median_L_inst_m": float(subset["L_inst"].median()),
                "median_obs_u10_mps": float(subset["obs_u10"].median()),
                "median_window_mean_u10_mps": float(subset["window_mean_u10"].median()),
                "median_pattern_lifetime_h": float(subset["pattern_lifetime_h"].median()),
                "median_n_events": float(subset["n_events"].median()),
                "median_comparison_bias_m": float(subset["comparison_bias_m"].median()),
                "mae_comparison_m": float(subset["comparison_bias_m"].abs().mean()),
            }
        )
    return pd.DataFrame(rows)


def plot_predicted_vs_observed(frame: pd.DataFrame, output_path: Path) -> None:
    """Scatter predicted wiggle spacing against observed wiggle spacing by site."""
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), sharex=True, sharey=True)
    layer_panels = [
        ("L_inst", "Onset width"),
        ("mechanical_capped_width", "Mechanical capped width"),
        ("comparison_spacing", "Comparison spacing"),
    ]
    for ax, (column, title) in zip(axes, layer_panels, strict=True):
        combined_min = min(frame["wiggle_spacing"].min(), frame[column].min())
        combined_max = max(frame["wiggle_spacing"].max(), frame[column].max())
        line = np.array([combined_min * 0.9, combined_max * 1.05])
        ax.plot(line, line, linestyle="--", linewidth=1.0, color="#666666")
        for site_id in SITE_ORDER:
            subset = frame.loc[frame["site_id"] == site_id]
            if subset.empty:
                continue
            ax.scatter(
                subset["wiggle_spacing"],
                subset[column],
                s=46,
                alpha=0.9,
                color=SITE_COLORS[site_id],
                edgecolors="white",
                linewidths=0.5,
                label=SITE_LABELS[site_id],
            )
        rho, _ = _spearman(frame, "wiggle_spacing", column)
        ax.set_title(f"{title}\nall-site rho={rho:.3f}" if math.isfinite(rho) else title)
        ax.set_xlabel("Observed wiggle spacing [m]")
        ax.grid(alpha=0.2, linewidth=0.6)
    axes[0].set_ylabel("Model spacing [m]")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=SITE_COLORS[site_id],
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=8,
            label=SITE_LABELS[site_id],
        )
        for site_id in SITE_ORDER
        if not frame.loc[frame["site_id"] == site_id].empty
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=max(1, len(handles)),
        frameon=False,
    )
    fig.suptitle("Cross-lake wiggle observations vs current model layers", y=0.98)
    fig.tight_layout(rect=(0, 0.08, 1, 0.92))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_temporal_development(frame: pd.DataFrame, output_path: Path) -> None:
    """Plot wiggle spacing development through time with matched wind diagnostics."""
    sites = [site_id for site_id in SITE_ORDER if not frame.loc[frame["site_id"] == site_id].empty]
    fig, axes = plt.subplots(len(sites), 1, figsize=(13.5, 3.4 * len(sites)), sharex=False)
    if len(sites) == 1:
        axes = [axes]
    for ax, site_id in zip(axes, sites, strict=True):
        subset = frame.loc[frame["site_id"] == site_id].sort_values("image_date").copy()
        dates = pd.to_datetime(subset["image_date"])
        ax.plot(
            dates,
            subset["wiggle_spacing"],
            color="#111111",
            linewidth=1.6,
            marker="o",
            markersize=5,
            label="Observed wiggle",
        )
        ax.plot(
            dates,
            subset["comparison_spacing"],
            color=SITE_COLORS[site_id],
            linewidth=1.4,
            marker="s",
            markersize=4,
            label="Comparison spacing",
        )
        ax.plot(
            dates,
            subset["L_inst"],
            color=SITE_COLORS[site_id],
            linewidth=1.1,
            linestyle="--",
            alpha=0.8,
            label="Onset width L_inst",
        )
        ax.set_ylabel("Spacing [m]")
        ax.set_title(SITE_LABELS[site_id])
        ax.grid(alpha=0.2, linewidth=0.6)
        ax2 = ax.twinx()
        ax2.bar(
            dates,
            subset["window_mean_u10"],
            width=20 if site_id == "taihu" else 12,
            alpha=0.18,
            color="#666666",
            label="Window mean U10",
        )
        ax2.plot(
            dates,
            subset["obs_u10"],
            color="#9a1f40",
            linewidth=1.2,
            marker="^",
            markersize=4,
            label="Observation U10",
        )
        ax2.set_ylabel("Wind speed [m/s]")
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
        if len(subset) <= 4:
            for _, row in subset.iterrows():
                ax.annotate(
                    row["case_id"],
                    (pd.Timestamp(row["image_date"]), float(row["wiggle_spacing"])),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=8,
                )
    primary_handles = [
        Line2D([0], [0], color="#111111", marker="o", linewidth=1.6, markersize=5, label="Observed wiggle"),
        Line2D([0], [0], color=SITE_COLORS[sites[0]], marker="s", linewidth=1.4, markersize=4, label="Comparison spacing"),
        Line2D([0], [0], color=SITE_COLORS[sites[0]], linestyle="--", linewidth=1.1, label="Onset width L_inst"),
        Line2D([0], [0], color="#9a1f40", marker="^", linewidth=1.2, markersize=4, label="Observation U10"),
    ]
    fig.legend(
        handles=primary_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.005),
        ncol=4,
        frameon=False,
    )
    fig.suptitle("Wiggle temporal development with matched wind diagnostics", y=0.985)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_wind_conditions(frame: pd.DataFrame, output_path: Path) -> None:
    """Plot observed wiggle spacing against wind metrics."""
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), sharey=True)
    specs = [
        ("obs_u10", "Observation U10 [m/s]"),
        ("window_mean_u10", "Window-mean U10 [m/s]"),
    ]
    for ax, (x_col, x_label) in zip(axes, specs, strict=True):
        for site_id in SITE_ORDER:
            subset = frame.loc[frame["site_id"] == site_id]
            if subset.empty:
                continue
            ax.scatter(
                subset[x_col],
                subset["wiggle_spacing"],
                s=40 + 18 * subset["n_events"],
                alpha=0.9,
                color=SITE_COLORS[site_id],
                edgecolors="white",
                linewidths=0.5,
                label=SITE_LABELS[site_id],
            )
        rho, _ = _spearman(frame, "wiggle_spacing", x_col)
        ax.set_title(f"all-site rho={rho:.3f}" if math.isfinite(rho) else "all-site rho=nan")
        ax.set_xlabel(x_label)
        ax.grid(alpha=0.2, linewidth=0.6)
    axes[0].set_ylabel("Observed wiggle spacing [m]")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=SITE_COLORS[site_id],
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=8,
            label=SITE_LABELS[site_id],
        )
        for site_id in SITE_ORDER
        if not frame.loc[frame["site_id"] == site_id].empty
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=max(1, len(handles)),
        frameon=False,
    )
    fig.suptitle("Wiggle spacing vs matched wind conditions", y=0.98)
    fig.tight_layout(rect=(0, 0.08, 1, 0.92))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_selected_hourly_evolution(
    frame: pd.DataFrame,
    observations: pd.DataFrame,
    cache_dir: Path,
    lookback_hours: float,
    output_path: Path,
) -> None:
    """Plot hourly wind evolution for Neagh wiggles and the three largest non-Neagh wiggles."""
    selected_ids = frame.loc[frame["site_id"] == "neagh", "case_id"].tolist()
    selected_ids += (
        frame.loc[frame["site_id"] != "neagh"]
        .sort_values("wiggle_spacing", ascending=False)
        .head(3)["case_id"]
        .tolist()
    )
    selected_ids = list(dict.fromkeys(selected_ids))
    selected = frame.set_index("case_id").loc[selected_ids].reset_index()
    obs_lookup = observations.set_index("observation_id", drop=False)

    fig, axes = plt.subplots(len(selected), 1, figsize=(12.8, 2.6 * len(selected)), sharex=True)
    if len(selected) == 1:
        axes = [axes]

    for ax, (_, row) in zip(axes, selected.iterrows(), strict=True):
        observation_row = obs_lookup.loc[row["case_id"]]
        matched = build_matched_case(
            observation_row=observation_row,
            cache_dir=cache_dir,
            site_registry=default_site_proxy_registry(),
            lookback_hours=lookback_hours,
        )
        weather = matched.weather_data.copy().sort_values("timestamp").reset_index(drop=True)
        rel_hours = (
            pd.to_datetime(weather["timestamp"], utc=True) - pd.Timestamp(matched.observation_time_utc)
        ).dt.total_seconds() / 3600.0
        color = SITE_COLORS[str(row["site_id"])]
        ax.plot(
            rel_hours,
            weather["U10"],
            color=color,
            linewidth=1.6,
            marker="o",
            markersize=4,
        )
        ax.scatter(
            [rel_hours.iloc[-1]],
            [weather["U10"].iloc[-1]],
            color="#9a1f40",
            s=40,
            zorder=4,
        )
        ax.axvline(0.0, color="#555555", linestyle="--", linewidth=1.0)
        ax.axhline(4.0, color="#888888", linestyle=":", linewidth=0.9)
        ax.set_ylabel("U10 [m/s]")
        ax.grid(alpha=0.2, linewidth=0.6)
        ax.set_title(
            f"{SITE_LABELS[str(row['site_id'])]}  {pd.Timestamp(row['image_date']).date()}  "
            f"obs={row['wiggle_spacing']:.1f} m  pred={row['comparison_spacing']:.1f} m"
        )
        ax.text(
            0.01,
            0.96,
            (
                f"window mean={row['window_mean_u10']:.2f} m/s, "
                f"obs U10={row['obs_u10']:.2f} m/s, "
                f"n_events={int(row['n_events'])}, "
                f"lifetime={row['pattern_lifetime_h']:.1f} h"
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 2.0},
        )

    axes[-1].set_xlabel("Hours to observation")
    handles = [
        Line2D([0], [0], color="#333333", linewidth=1.6, marker="o", markersize=4, label="Hourly U10"),
        Line2D([0], [0], color="#9a1f40", marker="o", linewidth=0.0, markersize=6, label="Observation-time U10"),
        Line2D([0], [0], color="#555555", linestyle="--", linewidth=1.0, label="Observation time"),
        Line2D([0], [0], color="#888888", linestyle=":", linewidth=0.9, label="4 m/s guide"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.005),
        ncol=4,
        frameon=False,
    )
    fig.suptitle(
        "Hourly wind evolution for Neagh wiggles and the three largest non-Neagh wiggles",
        y=0.99,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def print_summary(frame: pd.DataFrame, layer_summary: pd.DataFrame, site_summary: pd.DataFrame) -> None:
    """Print high-signal wiggle-only cross-lake summaries."""
    print(f"Wiggle-bearing observations analysed: {len(frame)}")
    print(
        "Site counts: "
        + ", ".join(
            f"{site_id}={int((frame['site_id'] == site_id).sum())}"
            for site_id in SITE_ORDER
            if (frame["site_id"] == site_id).any()
        )
    )
    overall = layer_summary.loc[
        (layer_summary["site_id"] == "all")
        & (layer_summary["x"].isin(["L_inst", "mechanical_capped_width", "comparison_spacing"]))
    ].copy()
    print("\nOverall wiggle layer summary")
    print(overall.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print("\nPer-site comparison-spacing summary")
    print(site_summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--comparison-output",
        type=Path,
        default=REPO_ROOT / "outputs" / "comparison_matched_current_code_rerun",
        help="Current matched comparison output directory [path]",
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "observations.csv",
        help="Benchmark observation CSV [path]",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "era5_cache",
        help="Matched ERA5/Open-Meteo cache directory [path]",
    )
    parser.add_argument(
        "--lookback-hours",
        type=float,
        default=6.0,
        help="Wind-history lookback used for matched cases [h]",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs" / "rank_audit" / "wiggle_crosslake_audit",
        help="Directory for wiggle-only audit outputs [path]",
    )
    return parser.parse_args()


def main() -> None:
    """Run the cross-lake wiggle audit and write tables/plots."""
    args = parse_args()
    case_dir = _resolve_case_dir(args.comparison_output)
    observations = load_wiggle_observations(args.observations)
    frame = build_wiggle_audit_table(
        case_dir=case_dir,
        observation_path=args.observations,
        cache_dir=args.cache_dir,
        lookback_hours=args.lookback_hours,
    )
    layer_summary = build_layer_summary(frame)
    site_summary = build_site_summary(frame)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "wiggle_crosslake_rank_audit.csv", index=False)
    layer_summary.to_csv(args.output_dir / "wiggle_crosslake_layer_summary.csv", index=False)
    site_summary.to_csv(args.output_dir / "wiggle_crosslake_site_summary.csv", index=False)

    plot_predicted_vs_observed(
        frame=frame,
        output_path=args.output_dir / "wiggle_predicted_vs_observed_by_site.png",
    )
    plot_temporal_development(
        frame=frame,
        output_path=args.output_dir / "wiggle_temporal_development_by_site.png",
    )
    plot_wind_conditions(
        frame=frame,
        output_path=args.output_dir / "wiggle_spacing_vs_wind_conditions.png",
    )
    plot_selected_hourly_evolution(
        frame=frame,
        observations=observations,
        cache_dir=args.cache_dir,
        lookback_hours=args.lookback_hours,
        output_path=args.output_dir / "selected_wiggle_hourly_wind_evolution.png",
    )

    print(f"Wrote wiggle cross-lake audit to {args.output_dir}")
    print_summary(frame=frame, layer_summary=layer_summary, site_summary=site_summary)


if __name__ == "__main__":
    main()
