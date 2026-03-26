#!/usr/bin/env python3
"""
Fetch-convention sensitivity analysis for the Neagh directional rerun.

This script reruns the Neagh-only analysis for both fetch conventions:
- `from`: use the meteorological wind-direction bearing directly for fetch lookup
- `to`:   use the downwind bearing (wind direction + 180 degrees)

Fetch bins are point-relative azimuths from the matched Neagh point. Example:
`fetch_0` means fetch from the point toward geographic north, so it is the
relevant bin for a northerly wind.

It uses the retained annotation `category` from the rebuilt point summary so
matched points remain distinguishable as:
- manual
- stream
- wiggle

Outputs:
- combined case table for both conventions
- per-convention and per-category summary metrics
- predicted-vs-observed sensitivity plots
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import warnings

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from neagh_directional_rerun import (  # type: ignore
    load_neagh_observations,
    match_points,
    rerun_observation,
)


CATEGORY_ORDER = ["manual", "stream", "wiggle"]
CATEGORY_COLORS = {
    "manual": "#1f5aa6",
    "stream": "#d1495b",
    "wiggle": "#e9c46a",
}
CONVENTION_COLORS = {
    "from": "#1f5aa6",
    "to": "#c65a1e",
}
CONVENTION_MARKERS = {
    "from": "o",
    "to": "^",
}

def _save_figure(fig: plt.Figure, output_path: Path) -> Path:
    """Save a Matplotlib figure and close it."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output_path


def rerun_both_conventions(
    observations_path: Path,
    point_summary_path: Path,
    era5_cache_dir: Path,
    lookback_hours: float,
) -> pd.DataFrame:
    """
    Rerun Neagh for both fetch conventions and return the combined case table.

    Parameters:
        observations_path: Base observation table [path]
        point_summary_path: Neagh point summary [path]
        era5_cache_dir: ERA5/Open-Meteo cache directory [path]
        lookback_hours: History window used for coarsening/disruption [h]

    Returns:
        Combined `from` + `to` rerun table [mixed]
    """
    observations = load_neagh_observations(observations_path)
    matched = match_points(observations, point_summary_path)
    rows: list[dict] = []
    for angle_mode in ("from", "to"):
        for _, row in matched.iterrows():
            record = rerun_observation(
                observation_row=row,
                cache_dir=era5_cache_dir,
                lookback_hours=lookback_hours,
                angle_mode=angle_mode,
                drag_method="coare35",
                drift_method="webb_fox_kemper",
            )
            record["fetch_angle_mode"] = angle_mode
            rows.append(record)
    combined = pd.DataFrame(rows).sort_values(["fetch_angle_mode", "case_id"]).reset_index(drop=True)
    combined["annotation_category"] = pd.Categorical(
        combined["annotation_category"],
        categories=CATEGORY_ORDER,
        ordered=True,
    )
    return combined


def build_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Build per-convention and per-category summary metrics.

    Parameters:
        frame: Combined sensitivity dataframe [mixed]

    Returns:
        Summary dataframe with rank and error metrics [mixed]
    """
    rows: list[dict] = []
    for angle_mode in ("from", "to"):
        for annotation_category in [None] + CATEGORY_ORDER:
            subset = frame.loc[frame["fetch_angle_mode"] == angle_mode].copy()
            if annotation_category is not None:
                subset = subset.loc[
                    subset["annotation_category"] == annotation_category
                ].copy()
            if subset.empty:
                continue
            observed = subset["observed_spacing"].to_numpy(dtype=float)
            predicted = subset["comparison_spacing"].to_numpy(dtype=float)
            if len(subset) >= 2 and np.unique(observed).size > 1 and np.unique(predicted).size > 1:
                rho, p_value = spearmanr(observed, predicted)
            else:
                rho = float("nan")
                p_value = float("nan")
            error = predicted - observed
            rows.append(
                {
                    "fetch_angle_mode": angle_mode,
                    "annotation_category": (
                        "all" if annotation_category is None else annotation_category
                    ),
                    "n": int(len(subset)),
                    "spearman_rho": float(rho),
                    "spearman_p_value": float(p_value),
                    "bias_m": float(np.mean(error)),
                    "mae_m": float(np.mean(np.abs(error))),
                    "median_observed_m": float(np.median(observed)),
                    "median_predicted_m": float(np.median(predicted)),
                    "median_point_mean_spacing_m": float(
                        np.median(subset["point_mean_spacing_m"])
                    ),
                }
            )
    return pd.DataFrame(rows)


def plot_predicted_vs_observed(frame: pd.DataFrame, output_path: Path) -> Path:
    """Write side-by-side predicted-vs-observed scatters for both conventions."""
    lower = float(min(frame["observed_spacing"].min(), frame["comparison_spacing"].min()))
    upper = float(max(frame["observed_spacing"].max(), frame["comparison_spacing"].max()))

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharex=True, sharey=True)
    for ax, angle_mode in zip(axes, ("from", "to"), strict=True):
        subset = frame.loc[frame["fetch_angle_mode"] == angle_mode].copy()
        observed = subset["observed_spacing"].to_numpy(dtype=float)
        predicted = subset["comparison_spacing"].to_numpy(dtype=float)
        rho, _ = spearmanr(observed, predicted) if len(subset) >= 2 else (float("nan"), float("nan"))
        for annotation_category in CATEGORY_ORDER:
            class_rows = subset.loc[
                subset["annotation_category"] == annotation_category
            ].copy()
            if class_rows.empty:
                continue
            ax.scatter(
                class_rows["observed_spacing"],
                class_rows["comparison_spacing"],
                s=56,
                color=CATEGORY_COLORS[annotation_category],
                marker=CONVENTION_MARKERS[angle_mode],
                alpha=0.85,
                edgecolors="white",
                linewidths=0.6,
            )
        ax.plot([lower, upper], [lower, upper], linestyle="--", color="#222222", linewidth=1.0)
        ax.set_title(f"Fetch = {angle_mode}  |  rho = {rho:.3f}")
        ax.set_xlabel("Observed spacing [m]")
        ax.grid(alpha=0.25, linewidth=0.6)
    axes[0].set_ylabel("Predicted spacing [m]")
    category_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.6,
            markersize=8,
            label=annotation_category,
        )
        for annotation_category, color in CATEGORY_COLORS.items()
    ]
    legend_category = fig.legend(
        handles=category_handles,
        loc="upper center",
        bbox_to_anchor=(0.50, 1.02),
        ncol=3,
        frameon=False,
        title="Annotation category",
    )
    fig.add_artist(legend_category)
    return _save_figure(fig, output_path)


def plot_error_by_stream_class(frame: pd.DataFrame, output_path: Path) -> Path:
    """Write prediction-error scatter/median plot grouped by annotation category."""
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    class_positions = {category: idx for idx, category in enumerate(CATEGORY_ORDER)}
    offsets = {"from": -0.16, "to": 0.16}

    rng = np.random.default_rng(20260325)
    for angle_mode in ("from", "to"):
        subset = frame.loc[frame["fetch_angle_mode"] == angle_mode].copy()
        subset["error_m"] = subset["comparison_spacing"] - subset["observed_spacing"]
        for annotation_category in CATEGORY_ORDER:
            x_center = class_positions[annotation_category] + offsets[angle_mode]
            all_class_rows = subset.loc[
                subset["annotation_category"] == annotation_category
            ].copy()
            if all_class_rows.empty:
                continue
            jitter = rng.uniform(-0.05, 0.05, size=len(all_class_rows))
            ax.scatter(
                np.full(len(all_class_rows), x_center) + jitter,
                all_class_rows["error_m"],
                s=42,
                color=CATEGORY_COLORS[annotation_category],
                marker=CONVENTION_MARKERS[angle_mode],
                alpha=0.8,
                edgecolors="white",
                linewidths=0.5,
            )
            median_error = float(all_class_rows["error_m"].median())
            ax.plot(
                [x_center - 0.08, x_center + 0.08],
                [median_error, median_error],
                color=CONVENTION_COLORS[angle_mode],
                linewidth=2.0,
            )
    ax.axhline(0.0, color="#444444", linestyle="--", linewidth=1.0)
    ax.set_xticks(range(len(CATEGORY_ORDER)), [label.title() for label in CATEGORY_ORDER])
    ax.set_ylabel("Predicted - observed [m]")
    ax.set_xlabel("Annotation category")
    ax.set_title("Prediction error by annotation category and fetch convention")
    convention_handles = [
        plt.Line2D(
            [0],
            [0],
            marker=CONVENTION_MARKERS[mode],
            color=CONVENTION_COLORS[mode],
            linestyle="",
            markersize=8,
            label=mode,
        )
        for mode in ("from", "to")
    ]
    category_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.6,
            markersize=8,
            label=annotation_category,
        )
        for annotation_category, color in CATEGORY_COLORS.items()
    ]
    legend_convention = ax.legend(
        handles=convention_handles,
        frameon=False,
        title="Fetch convention",
        loc="upper left",
    )
    ax.add_artist(legend_convention)
    ax.legend(
        handles=category_handles,
        frameon=False,
        title="Annotation category",
        loc="upper right",
    )
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    return _save_figure(fig, output_path)


def plot_prediction_shift(frame: pd.DataFrame, output_path: Path) -> Path:
    """Write the per-case prediction shift induced by changing fetch convention."""
    pivot = (
        frame.pivot(
            index="case_id",
            columns="fetch_angle_mode",
            values="comparison_spacing",
        )
        .rename_axis(columns=None)
        .reset_index()
    )
    base = frame.drop_duplicates("case_id")[
        ["case_id", "observed_spacing", "annotation_category", "point_mean_spacing_m"]
    ]
    merged = base.merge(pivot, on="case_id", how="left", validate="one_to_one")
    merged["delta_to_minus_from_m"] = merged["to"] - merged["from"]

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for annotation_category in CATEGORY_ORDER:
        class_rows = merged.loc[
            merged["annotation_category"] == annotation_category
        ].copy()
        if class_rows.empty:
            continue
        ax.scatter(
            class_rows["observed_spacing"],
            class_rows["delta_to_minus_from_m"],
            s=54,
            color=CATEGORY_COLORS[annotation_category],
            alpha=0.85,
            edgecolors="white",
            linewidths=0.6,
            label=annotation_category,
        )
    ax.axhline(0.0, color="#444444", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Observed spacing [m]")
    ax.set_ylabel("Predicted spacing shift: to - from [m]")
    ax.set_title("Fetch-convention sensitivity by observed spacing")
    category_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.6,
            markersize=8,
            label=annotation_category,
        )
        for annotation_category, color in CATEGORY_COLORS.items()
    ]
    ax.legend(handles=category_handles, frameon=False, title="Annotation category", loc="upper left")
    ax.grid(alpha=0.25, linewidth=0.6)
    return _save_figure(fig, output_path)


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
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs" / "rank_audit" / "neagh_fetch_convention_sensitivity",
        help="Directory for combined tables and plots [path]",
    )
    return parser.parse_args()


def main() -> None:
    """Run the fetch-convention sensitivity analysis and write tables/plots."""
    args = parse_args()
    warnings.filterwarnings(
        "ignore",
        message=r"Coarsened width .* exceeds cap .*",
        category=UserWarning,
    )

    combined = rerun_both_conventions(
        observations_path=args.observations,
        point_summary_path=args.point_summary,
        era5_cache_dir=args.era5_cache_dir,
        lookback_hours=args.lookback_hours,
    )
    summary = build_summary(combined)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = args.output_dir / "neagh_fetch_convention_sensitivity.csv"
    summary_path = args.output_dir / "neagh_fetch_convention_summary.csv"
    combined.to_csv(combined_path, index=False)
    summary.to_csv(summary_path, index=False)

    pred_obs_path = plot_predicted_vs_observed(
        combined,
        args.output_dir / "predicted_vs_observed_by_fetch_convention_and_category.png",
    )
    error_path = plot_error_by_stream_class(
        combined,
        args.output_dir / "prediction_error_by_annotation_category.png",
    )
    delta_path = plot_prediction_shift(
        combined,
        args.output_dir / "prediction_shift_to_minus_from_by_category.png",
    )

    print(f"Wrote combined table: {combined_path}")
    print(f"Wrote summary table: {summary_path}")
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print(f"Wrote plot: {pred_obs_path}")
    print(f"Wrote plot: {error_path}")
    print(f"Wrote plot: {delta_path}")


if __name__ == "__main__":
    main()
