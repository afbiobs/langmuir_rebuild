#!/usr/bin/env python3
"""
Layer-by-layer rank audit for the Langmuir rebuild.

KEY CONTEXT
- L_inst = 2*pi*depth / l_cNL  (raw onset width, pre-coarsening)
- Comparison spacing = L_inst * 2^n_events (after coarsening)
- Current full-set: Spearman rho = -0.464 (ANTI-correlated)
- Question: does L_inst already anti-correlate, or does coarsening flip it?
- 67 cases total, multiple lakes/sites
- ALSO include wiggle spacing from original observations in the analysis as the
  wiggle spacing may be capturing a different process

Implementation note:
- `coarsened_width` below is the raw post-coarsening hierarchy width from the
  CL diagnostics (`raw_coarsened_width_m`), before any mechanical cell-width
  cap. This isolates whether the sign flip happens in the coarsening layer or
  only after the final comparison/visibility mapping.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
from scipy.stats import spearmanr


AUDIT_COLUMNS = [
    "case_id",
    "observed_spacing",
    "U10",
    "depth",
    "Ra",
    "L_inst",
    "coarsened_width",
    "comparison_spacing",
    "pattern_lifetime_s",
    "u_star_water",
    "nu_T_vertical",
    "A_H",
    "site_id",
    "wiggle_spacing",
    "manual_spacing",
    "measurement_method",
    "n_events",
    "visible_n_events",
    "mechanical_capped_width",
    "cap_binding",
]

RANK_INPUT_COLUMNS = [
    "U10",
    "depth",
    "Ra",
    "L_inst",
    "coarsened_width",
    "comparison_spacing",
    "pattern_lifetime_s",
    "u_star_water",
    "nu_T_vertical",
    "A_H",
]

LAYER_CHAIN = ["U10", "Ra", "L_inst", "coarsened_width", "comparison_spacing"]


def _as_float(value: Any) -> float:
    """Convert scalar JSON values to float, preserving NaN for missing entries."""
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


def _finite(value: Any) -> bool:
    """Return True when `value` can be interpreted as a finite float."""
    return math.isfinite(_as_float(value))


def _first_finite(*values: Any) -> float:
    """Return the first finite numeric value from `values`, else NaN."""
    for value in values:
        number = _as_float(value)
        if math.isfinite(number):
            return number
    return float("nan")


def _as_int(value: Any, default: int = 0) -> int:
    """Convert a numeric value to int, returning `default` for missing values."""
    number = _as_float(value)
    if math.isfinite(number):
        return int(number)
    return int(default)


def _as_bool(value: Any, default: bool = False) -> bool:
    """Convert scalar JSON values to bool without treating NaN as True."""
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
    """Resolve either a comparison output dir or a direct case_diagnostics dir."""
    if path.name == "case_diagnostics":
        case_dir = path
    else:
        case_dir = path / "case_diagnostics"
    if not case_dir.exists():
        raise FileNotFoundError(f"Case diagnostics directory not found: {case_dir}")
    return case_dir


def load_observations(path: Path) -> pd.DataFrame:
    """
    Load the raw observation table and reconstruct the comparison target spacing.

    Parameters:
        path: Observation CSV path [path]

    Returns:
        Observation table indexed by observation_id with target/manual/wiggle
        spacings preserved [mixed]
    """
    frame = pd.read_csv(path)
    frame["target_spacing_m"] = frame["manual_spacing_m"].combine_first(
        frame["wiggle_spacing_m"]
    )
    return frame.set_index("observation_id", drop=False)


def _compute_l_inst(depth_m: Any, critical_result: dict[str, Any], base_width_m: Any) -> float:
    """
    Compute the onset width L_inst [m] from l_c when available.

    Parameters:
        depth_m:          Water depth h [m]
        critical_result:  CL critical solution diagnostics [-]
        base_width_m:     Recorded onset width from diagnostics [m]

    Returns:
        L_inst: Raw onset width before coarsening [m]
    """
    depth = _as_float(depth_m)
    l_c = _first_finite(
        critical_result.get("l_c"),
        critical_result.get("lcNL"),
    )
    if math.isfinite(depth) and depth > 0.0 and math.isfinite(l_c) and l_c > 0.0:
        return float((2.0 * math.pi * depth) / l_c)
    return _as_float(base_width_m)


def extract_case_row(
    case_file: Path,
    observations: pd.DataFrame,
    model: str,
) -> dict[str, Any]:
    """
    Extract one CL case into a flat audit row.

    Parameters:
        case_file:     Per-case diagnostic JSON [path]
        observations:  Observation lookup table [mixed]
        model:         Candidate model key within the case JSON [-]

    Returns:
        Flat audit row suitable for CSV export and rank analysis [mixed]
    """
    with case_file.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)

    model_payload = payload["models"][model]
    prediction_row = model_payload.get("prediction_row", {})
    full = model_payload.get("full_diagnostic", {})
    case_input = payload.get("case_input", {})
    forcing = full.get("forcing_summary", {})
    coarsening = full.get("coarsening", {})
    intermediate = full.get("intermediate", {})
    critical_result = intermediate.get("critical_result", {})

    case_id = str(
        payload.get("observation_id")
        or prediction_row.get("observation_id")
        or case_input.get("observation_id")
    )
    if case_id not in observations.index:
        raise KeyError(f"{case_id} is missing from the observation table.")
    observation = observations.loc[case_id]

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
    visible_n_events = _as_int(
        _first_finite(coarsening.get("visible_n_events"), prediction_row.get("n_visible_events"))
    )
    coarsened_width = _first_finite(
        coarsening.get("raw_coarsened_width_m"),
        onset_width * (2.0 ** n_events) if math.isfinite(onset_width) else float("nan"),
    )
    comparison_spacing = _first_finite(
        coarsening.get("visible_spacing_m"),
        prediction_row.get("comparison_spacing_m"),
        prediction_row.get("raw_predicted_spacing_m"),
    )
    pattern_lifetime = _first_finite(
        full.get("pattern_lifetime_s"),
        case_input.get("pattern_lifetime_s"),
        prediction_row.get("pattern_lifetime_s"),
    )
    return {
        "case_id": case_id,
        "observed_spacing": _first_finite(
            observation.get("target_spacing_m"),
            case_input.get("observed_spacing_m"),
            prediction_row.get("observed_spacing_m"),
        ),
        "U10": _first_finite(
            forcing.get("U10"),
            case_input.get("representative_u10_mps"),
            prediction_row.get("representative_u10_mps"),
        ),
        "depth": depth,
        "Ra": _first_finite(full.get("Ra"), prediction_row.get("Ra")),
        "L_inst": _compute_l_inst(depth, critical_result, onset_width),
        "coarsened_width": coarsened_width,
        "comparison_spacing": comparison_spacing,
        "pattern_lifetime_s": pattern_lifetime,
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
        "site_id": str(
            case_input.get("site_id")
            or prediction_row.get("site_id")
            or "unknown"
        ),
        "wiggle_spacing": _as_float(observation.get("wiggle_spacing_m")),
        "manual_spacing": _as_float(observation.get("manual_spacing_m")),
        "measurement_method": str(
            case_input.get("measurement_method")
            or prediction_row.get("measurement_method")
            or (
                "manual"
                if _finite(observation.get("manual_spacing_m"))
                else "wiggle"
            )
        ),
        "n_events": n_events,
        "visible_n_events": visible_n_events,
        "mechanical_capped_width": _first_finite(coarsening.get("coarsened_width_m")),
        "cap_binding": _as_bool(
            coarsening.get("cap_binding", prediction_row.get("cap_binding", False))
        ),
    }


def build_rank_audit_table(
    case_dir: Path,
    observation_path: Path,
    model: str,
) -> pd.DataFrame:
    """
    Build the flat rank-audit table from case diagnostics and observations.

    Parameters:
        case_dir:           Directory containing `case_*.json` files [path]
        observation_path:   Observation CSV path [path]
        model:              Candidate model to audit [-]

    Returns:
        Rank-audit table with one row per case [mixed]
    """
    observations = load_observations(observation_path)
    rows = [
        extract_case_row(case_file, observations=observations, model=model)
        for case_file in sorted(case_dir.glob("case_*.json"))
    ]
    if not rows:
        raise FileNotFoundError(f"No case JSON files found under {case_dir}")

    frame = pd.DataFrame(rows)
    frame = frame[AUDIT_COLUMNS].sort_values("case_id").reset_index(drop=True)
    return frame


def _spearman_summary(frame: pd.DataFrame, target_col: str, x_cols: list[str]) -> pd.DataFrame:
    """
    Compute Spearman rank correlation for one target against multiple inputs.

    Parameters:
        frame:       Audit dataframe [mixed]
        target_col:  Target spacing column [m]
        x_cols:      Candidate predictor columns [mixed]

    Returns:
        Summary table with n, rho, and p-value per predictor [mixed]
    """
    rows: list[dict[str, Any]] = []
    for column in x_cols:
        valid = frame[[target_col, column]].dropna()
        if len(valid) < 2 or valid[target_col].nunique() < 2 or valid[column].nunique() < 2:
            rho = float("nan")
            p_value = float("nan")
        else:
            rho, p_value = spearmanr(valid[target_col], valid[column])
            rho = _as_float(rho)
            p_value = _as_float(p_value)
        rows.append(
            {
                "x": column,
                "n": int(len(valid)),
                "rho": rho,
                "p_value": p_value,
            }
        )
    return pd.DataFrame(rows)


def _first_negative_layer(summary: pd.DataFrame, chain: list[str]) -> str:
    """Return the earliest column in `chain` with a negative finite Spearman rho."""
    indexed = summary.set_index("x")
    for column in chain:
        if column not in indexed.index:
            continue
        rho = _as_float(indexed.at[column, "rho"])
        if math.isfinite(rho) and rho < 0.0:
            return f"{column} (rho={rho:.3f})"
    return "none"


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
            "Layer flip: hydrodynamic control remains non-negative at Ra "
            f"(rho={ra_rho:.3f}), but onset width L_inst is already anti-correlated "
            f"(rho={onset_rho:.3f})."
        )
    if math.isfinite(onset_rho) and onset_rho < 0.0:
        return (
            "Layer flip: onset width L_inst is already anti-correlated "
            f"(rho={onset_rho:.3f}); coarsening changes it to "
            f"{coarsened_rho:.3f} and final comparison spacing to {comparison_rho:.3f}."
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
            "Layer flip: onset and coarsening stay non-negative, and the final "
            f"comparison layer flips the sign (rho={comparison_rho:.3f})."
        )
    return (
        "Layer flip: no negative correlation appears in the onset/coarsening/"
        "comparison chain."
    )


def _print_summary(title: str, summary: pd.DataFrame) -> None:
    """Print one Spearman summary block."""
    print(f"\n{title}")
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))


def print_rank_audit(frame: pd.DataFrame) -> None:
    """
    Print full-set and per-site rank summaries.

    Parameters:
        frame: Rank-audit dataframe [mixed]

    Returns:
        None
    """
    print(f"Cases loaded: {len(frame)}")
    print(f"Sites: {', '.join(sorted(frame['site_id'].dropna().unique()))}")

    target_summary = _spearman_summary(frame, "observed_spacing", RANK_INPUT_COLUMNS)
    _print_summary("Full set: target observed_spacing", target_summary)
    print(
        "First negative layer in chain "
        f"{LAYER_CHAIN}: {_first_negative_layer(target_summary, LAYER_CHAIN)}"
    )
    print(_layer_flip_message(target_summary))

    wiggle_frame = frame.dropna(subset=["wiggle_spacing"]).copy()
    wiggle_summary = _spearman_summary(wiggle_frame, "wiggle_spacing", RANK_INPUT_COLUMNS)
    _print_summary(
        "Full set: wiggle_spacing subset",
        wiggle_summary,
    )
    print(
        "First negative layer in chain for wiggle spacing "
        f"{LAYER_CHAIN}: {_first_negative_layer(wiggle_summary, LAYER_CHAIN)}"
    )
    print(_layer_flip_message(wiggle_summary))

    for site_id, site_frame in frame.groupby("site_id", sort=True):
        site_summary = _spearman_summary(site_frame, "observed_spacing", RANK_INPUT_COLUMNS)
        _print_summary(
            f"Per-site target observed_spacing: {site_id} (n={len(site_frame)})",
            site_summary,
        )
        print(
            "First negative layer in chain "
            f"{LAYER_CHAIN}: {_first_negative_layer(site_summary, LAYER_CHAIN)}"
        )
        print(_layer_flip_message(site_summary))

        site_wiggle_frame = site_frame.dropna(subset=["wiggle_spacing"]).copy()
        if not site_wiggle_frame.empty:
            site_wiggle_summary = _spearman_summary(
                site_wiggle_frame,
                "wiggle_spacing",
                RANK_INPUT_COLUMNS,
            )
            _print_summary(
                f"Per-site wiggle_spacing: {site_id} (n={len(site_wiggle_frame)})",
                site_wiggle_summary,
            )
            print(
                "First negative layer in chain "
                f"{LAYER_CHAIN}: {_first_negative_layer(site_wiggle_summary, LAYER_CHAIN)}"
            )
            print(_layer_flip_message(site_wiggle_summary))


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the rank-audit script."""
    repo_root = Path(__file__).resolve().parents[2]
    default_output_dir = repo_root / "outputs" / "comparison_matched_current_code_rerun"
    default_observations = repo_root / "data" / "raw" / "observations.csv"
    audit_output_dir = repo_root / "outputs" / "rank_audit"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--comparison-dir",
        type=Path,
        default=default_output_dir,
        help="Comparison output directory containing case_diagnostics/ [path]",
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=default_observations,
        help="Observation CSV used to recover manual and wiggle spacing [path]",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="cl",
        help="Model key within each case diagnostic JSON [-]",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=audit_output_dir / "rank_audit.csv",
        help="CSV path for the flattened audit table [path]",
    )
    return parser.parse_args()


def main() -> None:
    """Run the rank audit and write the flattened case table to CSV."""
    args = parse_args()
    case_dir = _resolve_case_dir(args.comparison_dir)
    audit = build_rank_audit_table(
        case_dir=case_dir,
        observation_path=args.observations,
        model=args.model,
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(audit)} rows to {args.output_csv}")
    print_rank_audit(audit)


if __name__ == "__main__":
    main()
