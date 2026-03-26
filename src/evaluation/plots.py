"""
WP-06 comparison plotting helpers.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from src.evaluation.metrics import attractor_test


MODEL_ORDER = [
    "cl",
    "scaling",
    "baseline_constant",
    "baseline_linear_wind",
    "baseline_depth_scaled",
]

MODEL_LABELS = {
    "cl": "CL nonlinear",
    "scaling": "Scaling laws",
    "baseline_constant": "Baseline: constant",
    "baseline_linear_wind": "Baseline: wind",
    "baseline_depth_scaled": "Baseline: wind + depth",
}

MODEL_COLORS = {
    "cl": "#1f5aa6",
    "scaling": "#c65a1e",
    "baseline_constant": "#6d6d6d",
    "baseline_linear_wind": "#4d8f47",
    "baseline_depth_scaled": "#8c4fb6",
}


def _save_figure(fig: plt.Figure, output_path: str | Path) -> Path:
    """Save a Matplotlib figure to disk and close it."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return path


def _model_prediction_rows(prediction_table: pd.DataFrame, model: str) -> pd.DataFrame:
    """Return the comparison rows for one model."""
    rows = prediction_table.loc[prediction_table["model"] == model].copy()
    if rows.empty:
        raise ValueError(f"No comparison rows found for model '{model}'.")
    return rows.sort_values("observed_spacing_m").reset_index(drop=True)


def write_predicted_vs_observed_plot(
    prediction_table: pd.DataFrame,
    model: str,
    output_path: str | Path,
) -> Path:
    """
    Write a predicted-vs-observed spacing scatter plot.

    Parameters:
        prediction_table: Comparison rows with observed and predicted spacing [m]
        model:            Model identifier [-]
        output_path:      PNG destination path [-]

    Returns:
        Path to the saved plot [-]
    """
    rows = _model_prediction_rows(prediction_table, model)
    observed = rows["observed_spacing_m"].to_numpy(dtype=float)
    predicted = rows["comparison_spacing_m"].to_numpy(dtype=float)
    lower = float(min(np.min(observed), np.min(predicted)))
    upper = float(max(np.max(observed), np.max(predicted)))

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(
        observed,
        predicted,
        s=48,
        color=MODEL_COLORS.get(model, "#444444"),
        alpha=0.85,
        edgecolors="white",
        linewidths=0.5,
    )
    ax.plot([lower, upper], [lower, upper], linestyle="--", color="#222222", linewidth=1.0)
    ax.set_title(f"Predicted vs observed: {MODEL_LABELS.get(model, model)}")
    ax.set_xlabel("Observed spacing [m]")
    ax.set_ylabel("Predicted spacing [m]")
    ax.grid(alpha=0.25, linewidth=0.6)
    return _save_figure(fig, output_path)


def write_spacing_vs_wind_plot(
    prediction_table: pd.DataFrame,
    model: str,
    output_path: str | Path,
) -> Path:
    """
    Write a spacing-versus-wind diagnostic for one model.

    Parameters:
        prediction_table: Comparison rows with wind speed [m/s] and spacing [m]
        model:            Model identifier [-]
        output_path:      PNG destination path [-]

    Returns:
        Path to the saved plot [-]
    """
    rows = _model_prediction_rows(prediction_table, model).sort_values("representative_u10_mps")
    wind = rows["representative_u10_mps"].to_numpy(dtype=float)
    observed = rows["observed_spacing_m"].to_numpy(dtype=float)
    predicted = rows["comparison_spacing_m"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    for x_value, observed_value, predicted_value in zip(wind, observed, predicted):
        ax.plot(
            [x_value, x_value],
            [observed_value, predicted_value],
            color="#c9c9c9",
            linewidth=1.0,
            zorder=1,
        )
    ax.scatter(
        wind,
        observed,
        s=42,
        marker="o",
        facecolors="none",
        edgecolors="#1f1f1f",
        linewidths=1.0,
        label="Observed",
        zorder=3,
    )
    ax.scatter(
        wind,
        predicted,
        s=44,
        marker="^",
        color=MODEL_COLORS.get(model, "#444444"),
        alpha=0.9,
        label="Predicted",
        zorder=4,
    )
    ax.set_title(f"Spacing vs wind: {MODEL_LABELS.get(model, model)}")
    ax.set_xlabel("Representative U10 [m/s]")
    ax.set_ylabel("Spacing [m]")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25, linewidth=0.6)
    return _save_figure(fig, output_path)


def write_dynamic_range_comparison(
    metrics_table: pd.DataFrame,
    output_path: str | Path,
    subset_id: str = "BM-A_full",
) -> Path:
    """
    Write the full-set dynamic-range comparison plot.

    Parameters:
        metrics_table: Table with range coverage [-] and p90/p10 ratio [-]
        output_path:   PNG destination path [-]
        subset_id:     Benchmark subset identifier [-]

    Returns:
        Path to the saved plot [-]
    """
    rows = metrics_table.loc[metrics_table["subset_id"] == subset_id].copy()
    rows["model"] = pd.Categorical(rows["model"], MODEL_ORDER, ordered=True)
    rows = rows.sort_values("model")
    labels = [MODEL_LABELS.get(model, model) for model in rows["model"]]
    coverage = rows["range_coverage_fraction"].to_numpy(dtype=float)
    ratios = rows["dynamic_range_p90_p10"].to_numpy(dtype=float)
    y_pos = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    colors = [MODEL_COLORS.get(model, "#444444") for model in rows["model"]]
    ax.barh(y_pos, coverage, color=colors, alpha=0.85)
    ax.set_yticks(y_pos, labels)
    ax.set_xlabel("Observed-range coverage [-]")
    ax.set_title("Dynamic range comparison")
    ax.set_xlim(0.0, max(1.0, float(np.nanmax(coverage)) * 1.15))
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    for idx, (coverage_value, ratio_value) in enumerate(zip(coverage, ratios)):
        ratio_text = "nan" if np.isnan(ratio_value) else f"{ratio_value:.2f}"
        ax.text(
            coverage_value + 0.02,
            idx,
            f"p90/p10 = {ratio_text}",
            va="center",
            fontsize=8,
        )
    return _save_figure(fig, output_path)


def write_tail_coverage_comparison(
    metrics_table: pd.DataFrame,
    output_path: str | Path,
    subset_id: str = "BM-A_full",
) -> Path:
    """
    Write the full-set tail-coverage comparison plot.

    Parameters:
        metrics_table: Table with low/high/overall tail coverage [-]
        output_path:   PNG destination path [-]
        subset_id:     Benchmark subset identifier [-]

    Returns:
        Path to the saved plot [-]
    """
    rows = metrics_table.loc[metrics_table["subset_id"] == subset_id].copy()
    rows["model"] = pd.Categorical(rows["model"], MODEL_ORDER, ordered=True)
    rows = rows.sort_values("model")
    labels = [MODEL_LABELS.get(model, model) for model in rows["model"]]
    x_pos = np.arange(len(rows))
    width = 0.24

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.bar(
        x_pos - width,
        rows["tail_coverage_low"].fillna(0.0).to_numpy(dtype=float),
        width=width,
        color="#4d8f47",
        label="Low tail",
    )
    ax.bar(
        x_pos,
        rows["tail_coverage_high"].fillna(0.0).to_numpy(dtype=float),
        width=width,
        color="#c65a1e",
        label="High tail",
    )
    ax.bar(
        x_pos + width,
        rows["tail_coverage_overall"].fillna(0.0).to_numpy(dtype=float),
        width=width,
        color="#1f5aa6",
        label="Overall",
    )
    ax.set_xticks(x_pos, labels, rotation=20, ha="right")
    ax.set_ylabel("Tail coverage [-]")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Tail coverage comparison")
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False, ncol=3)
    return _save_figure(fig, output_path)


def _densest_band(
    predicted: np.ndarray,
    observed_range: tuple[float, float],
    band_fraction: float = 0.20,
) -> tuple[float, float, float]:
    """Locate the most concentrated prediction band over the observed range."""
    obs_min, obs_max = float(observed_range[0]), float(observed_range[1])
    band_width = band_fraction * (obs_max - obs_min)
    centres = np.linspace(obs_min + band_width / 2, obs_max - band_width / 2, 200)
    max_fraction = 0.0
    band_limits = (obs_min, obs_min + band_width)
    for centre in centres:
        lower = centre - band_width / 2
        upper = centre + band_width / 2
        fraction = float(np.mean((predicted >= lower) & (predicted <= upper)))
        if fraction >= max_fraction:
            max_fraction = fraction
            band_limits = (lower, upper)
    return band_limits[0], band_limits[1], max_fraction


def write_attractor_diagnostic(
    prediction_table: pd.DataFrame,
    model: str,
    output_path: str | Path,
) -> Path:
    """
    Write the attractor diagnostic for one model.

    Parameters:
        prediction_table: Comparison rows with predicted and observed spacing [m]
        model:            Model identifier [-]
        output_path:      PNG destination path [-]

    Returns:
        Path to the saved plot [-]
    """
    rows = _model_prediction_rows(prediction_table, model)
    predicted = rows["comparison_spacing_m"].to_numpy(dtype=float)
    observed = rows["observed_spacing_m"].to_numpy(dtype=float)
    observed_range = (float(np.min(observed)), float(np.max(observed)))
    attractor = attractor_test(predicted, observed_range=observed_range)
    band_lower, band_upper, max_fraction = _densest_band(predicted, observed_range)

    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.hist(
        predicted,
        bins=min(12, max(4, len(predicted))),
        color=MODEL_COLORS.get(model, "#444444"),
        alpha=0.8,
        edgecolor="white",
        linewidth=0.8,
    )
    ax.axvspan(band_lower, band_upper, color="#f0d8a8", alpha=0.35)
    ax.axvline(observed_range[0], color="#222222", linestyle="--", linewidth=1.0)
    ax.axvline(observed_range[1], color="#222222", linestyle="--", linewidth=1.0)
    ax.set_title(f"Attractor diagnostic: {MODEL_LABELS.get(model, model)}")
    ax.set_xlabel("Predicted spacing [m]")
    ax.set_ylabel("Case count [-]")
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    ax.text(
        0.02,
        0.95,
        (
            f"max fraction in 20% band = {max_fraction:.2%}\n"
            f"status = {'pass' if attractor['passed'] else 'fail'}"
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "#d7d7d7", "alpha": 0.9},
    )
    return _save_figure(fig, output_path)


def write_enhancement_index_timeseries(
    timeline_frame: pd.DataFrame,
    observation_id: str,
    output_path: str | Path,
) -> Path:
    """
    Write the enhancement-index time series for one representative case.

    Parameters:
        timeline_frame: Rows with timestamp [UTC] and development index [-]
        observation_id: Observation identifier [-]
        output_path:     PNG destination path [-]

    Returns:
        Path to the saved plot [-]
    """
    if timeline_frame.empty:
        raise ValueError("timeline_frame must contain at least one row.")

    frame = timeline_frame.copy()
    frame["timestamp"] = pd.to_datetime(frame["timestamp"], utc=True)
    frame = frame.sort_values("timestamp")

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for model in ("cl", "scaling"):
        rows = frame.loc[frame["model"] == model]
        if rows.empty:
            continue
        color = MODEL_COLORS.get(model, "#444444")
        if len(rows) == 1:
            ax.scatter(
                rows["timestamp"],
                rows["development_index"],
                s=48,
                color=color,
                label=MODEL_LABELS.get(model, model),
            )
        else:
            ax.plot(
                rows["timestamp"],
                rows["development_index"],
                marker="o",
                linewidth=1.6,
                markersize=3.8,
                color=color,
                label=MODEL_LABELS.get(model, model),
            )
    ax.set_title(f"Enhancement index timeseries: {observation_id}")
    ax.set_xlabel("Timestamp [UTC]")
    ax.set_ylabel("Development index [-]")
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False)
    fig.autofmt_xdate(rotation=25, ha="right")
    return _save_figure(fig, output_path)
