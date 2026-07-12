"""Plot station-by-station comparisons of SCHISM and downloaded USGS data.

The SCHISM file must contain elapsed model time in its first column and one
station per remaining column, in the same order as the accompanying BP file.
The observation directory must be an output from download_usgs_polygon.py.
"""

from __future__ import annotations

import argparse
import json
import re
from collections.abc import Mapping
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Extract.extract_schism_output import DEFAULT_EXTRACT_SCRIPT, ensure_model_file
from Stats.model_obs_io import (
    TIME_UNIT_TO_SECONDS,
    normalize_variable,
    read_bp_stations,
    read_model,
    read_observations,
)

MODEL_COLORS = ["tab:blue", "tab:green", "tab:orange"]
OBSERVATION_COLOR = "tab:red"
MAX_MODELS = len(MODEL_COLORS)
OVERVIEW_FIGSIZE_INCHES = (16, 9)
OVERVIEW_DPI = 240


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "station"


def _plot_station_axis(
    axis,
    station,
    models: Mapping[str, pd.DataFrame],
    observations: pd.DataFrame,
    variable: str,
    ylabel: str | None,
    start_time: pd.Timestamp | None,
    end_time: pd.Timestamp | None,
    *,
    overview: bool,
    show_legend: bool,
) -> None:
    station_id = str(station.station_id)
    station_obs = observations[observations["station_id"] == station_id]
    if start_time is not None:
        station_obs = station_obs[station_obs["time"] >= start_time]
    if end_time is not None:
        station_obs = station_obs[station_obs["time"] <= end_time]

    provenance_column = "derivation" if "derivation" in station_obs.columns else None
    groups = (
        station_obs.groupby(provenance_column, dropna=False)
        if provenance_column
        else [("observed", station_obs)]
    )
    for provenance, group in groups:
        label = "USGS" if overview else (
            f"USGS ({provenance})" if provenance else "USGS"
        )
        axis.plot(
            group["time"], group["value"], ".", color=OBSERVATION_COLOR,
            markersize=2 if overview else 3, label=label,
        )

    for (model_name, model), color in zip(models.items(), MODEL_COLORS):
        model_series = model[station_id]
        if start_time is not None:
            model_series = model_series[model_series.index >= start_time]
        if end_time is not None:
            model_series = model_series[model_series.index <= end_time]
        axis.plot(
            model_series.index, model_series, color=color,
            linewidth=0.8 if overview else 1.1, label=model_name,
        )

    station_name = str(station.station_name).strip()
    title = station_id if overview else (
        f"{station_id}: {station_name}" if station_name else station_id
    )
    axis.set_title(title, fontsize=9 if overview else None)
    if not overview:
        axis.set_xlabel("Time (UTC)")
    if ylabel:
        axis.set_ylabel(ylabel)
    else:
        units = station_obs.get("unit", pd.Series(dtype=str)).dropna().unique()
        unit_text = f" ({units[0]})" if len(units) == 1 else ""
        axis.set_ylabel(f"{variable}{unit_text}")
    axis.grid(alpha=0.3)
    if overview:
        axis.tick_params(axis="x", labelrotation=25, labelsize=7)
        axis.tick_params(axis="y", labelsize=7)
    else:
        axis.tick_params(axis="x", labelrotation=25)
    if show_legend:
        if overview:
            axis.legend(fontsize="x-small")
        else:
            axis.legend()


def plot_comparisons(
    models: Mapping[str, pd.DataFrame],
    observations: pd.DataFrame,
    stations: pd.DataFrame,
    output_dir: Path,
    variable: str,
    ylabel: str | None = None,
    start: str | None = None,
    end: str | None = None,
    dpi: int = 150,
    overview_layout: tuple[int, int] = (6, 3),
    overview_ylim: tuple[float, float] | None = None,
) -> list[Path]:
    """Write individual plots plus paginated overview plots."""

    if not models:
        raise ValueError("At least one model is required")
    if len(models) > MAX_MODELS:
        raise ValueError(
            f"At most {MAX_MODELS} models may be plotted; received {len(models)}"
        )
    overview_rows, overview_columns = overview_layout
    if overview_rows < 1 or overview_columns < 1:
        raise ValueError("Overview rows and columns must be positive")
    if overview_ylim is not None and overview_ylim[0] >= overview_ylim[1]:
        raise ValueError("Overview y-axis minimum must be less than its maximum")
    output_dir.mkdir(parents=True, exist_ok=True)
    start_time = pd.to_datetime(start, utc=True).tz_localize(None) if start else None
    end_time = pd.to_datetime(end, utc=True).tz_localize(None) if end else None
    written = []

    for station in stations.itertuples(index=False):
        station_id = str(station.station_id)
        fig, axis = plt.subplots(figsize=(10, 4))
        _plot_station_axis(
            axis, station, models, observations, variable, ylabel,
            start_time, end_time, overview=False, show_legend=True,
        )
        fig.autofmt_xdate()
        fig.tight_layout()
        path = output_dir / f"{safe_filename(station_id)}_{safe_filename(variable)}.png"
        fig.savefig(path, dpi=dpi)
        plt.close(fig)
        written.append(path)

    station_rows = list(stations.itertuples(index=False))
    page_size = overview_rows * overview_columns
    for page_number, offset in enumerate(
        range(0, len(station_rows), page_size), start=1
    ):
        page = station_rows[offset:offset + page_size]
        fig, axes = plt.subplots(
            overview_rows,
            overview_columns,
            figsize=OVERVIEW_FIGSIZE_INCHES,
            squeeze=False,
            sharex=True,
        )
        flat_axes = axes.ravel()
        for index, (axis, station) in enumerate(zip(flat_axes, page)):
            _plot_station_axis(
                axis, station, models, observations, variable, ylabel,
                start_time, end_time, overview=True, show_legend=index == 0,
            )
            if overview_ylim is not None:
                axis.set_ylim(*overview_ylim)
            if index % overview_columns != 0:
                axis.set_ylabel("")
        for axis in flat_axes[len(page):]:
            axis.set_visible(False)
        first = offset + 1
        last = offset + len(page)
        fig.suptitle(
            f"{variable}: USGS vs models, stations {first}-{last}", fontsize=15
        )
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        path = output_dir / (
            f"overview_{safe_filename(variable)}_{page_number:02d}.png"
        )
        fig.savefig(path, dpi=OVERVIEW_DPI)
        plt.close(fig)
        written.append(path)
    return written


def _load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def _comparison_items(config: dict, names: list[str] | None = None) -> list[dict]:
    if "comparisons" in config:
        comparisons = config["comparisons"]
    else:
        comparisons = [config]
    if not isinstance(comparisons, list) or not comparisons:
        raise ValueError("Config must contain a nonempty comparisons list")
    if names:
        requested = set(names)
        comparisons = [item for item in comparisons if item.get("name") in requested]
        missing = requested - {item.get("name") for item in comparisons}
        if missing:
            raise ValueError(f"Requested comparisons not found: {sorted(missing)}")
    return comparisons


def _path(value: str | Path) -> Path:
    return Path(value)


def _models_from_config(
    comparison: dict,
    extract_missing: bool,
    extract_mpi_np: int,
    extract_script: str | None,
    extract_config_dir: Path | None,
) -> tuple[list[Path], list[str]]:
    runs = comparison.get("runs")
    if not isinstance(runs, list) or not runs:
        raise ValueError(f"Comparison {comparison.get('name', '<unnamed>')} must define nonempty runs")
    labels = []
    files = []
    default_extract_config_dir = extract_config_dir or Path(comparison.get("output_dir", ".")) / "extract_configs"
    for run in runs:
        if not isinstance(run, dict):
            raise ValueError("Each run must be a dictionary")
        run_spec = {
            **run,
            "elev_out_file": run.get("elev_out_file") or run.get("model_file"),
            "station_bp_file": run.get("station_bp_file") or comparison["bp_file"],
        }
        if not run_spec.get("elev_out_file"):
            raise ValueError(f"Run {run.get('name', '<unnamed>')} is missing elev_out_file/model_file")
        ensure_model_file(
            run_spec,
            extract_missing=extract_missing,
            extract_mpi_np=int(run.get("extract_mpi_np", extract_mpi_np)),
            extract_script=extract_script or run.get("extract_script") or DEFAULT_EXTRACT_SCRIPT,
            extract_config_dir=run.get("extract_config_dir", default_extract_config_dir),
        )
        files.append(_path(run_spec["elev_out_file"]))
        labels.append(str(run.get("name") or run.get("label") or f"Model {len(labels) + 1}"))
    return files, labels


def run_comparison(
    obs_dir: Path,
    model_files: list[Path],
    model_labels: list[str] | None,
    bp_file: Path,
    model_start: str,
    variable: str,
    output_dir: Path,
    model_time_unit: str = "days",
    missing_value: float = -9999.0,
    ylabel: str | None = None,
    start: str | None = None,
    end: str | None = None,
    dpi: int = 150,
    overview_layout: tuple[int, int] = (6, 3),
    overview_ylim: tuple[float, float] | None = None,
) -> list[Path]:
    variable = normalize_variable(variable)
    stations = read_bp_stations(bp_file)
    observations = read_observations(obs_dir, variable)
    if len(model_files) > MAX_MODELS:
        raise ValueError(f"At most {MAX_MODELS} model files are supported")
    labels = model_labels or [
        f"Model {index}" for index in range(1, len(model_files) + 1)
    ]
    if len(labels) != len(model_files):
        raise ValueError("Provide exactly one model label for each model file, or omit all labels")
    if len(set(labels)) != len(labels):
        raise ValueError("Model labels must be unique")
    models = {
        label: read_model(
            model_file,
            stations,
            model_start,
            model_time_unit,
            missing_value,
        )
        for label, model_file in zip(labels, model_files)
    }
    written = plot_comparisons(
        models,
        observations,
        stations,
        output_dir,
        variable,
        ylabel=ylabel,
        start=start,
        end=end,
        dpi=dpi,
        overview_layout=tuple(overview_layout),
        overview_ylim=tuple(overview_ylim) if overview_ylim is not None else None,
    )
    first_model = next(iter(models.values()))
    matched = set(observations["station_id"]).intersection(first_model.columns)
    print(
        f"Wrote {len(written)} plots to {output_dir}; "
        f"{len(matched)} stations have observations."
    )
    return written


def run_config(
    config_file: Path,
    comparison_names: list[str] | None = None,
    extract_missing: bool = True,
    extract_mpi_np: int = 4,
    extract_script: str | None = None,
    extract_config_dir: Path | None = None,
) -> None:
    config = _load_json(config_file)
    for comparison in _comparison_items(config, comparison_names):
        model_files, model_labels = _models_from_config(
            comparison,
            extract_missing=extract_missing,
            extract_mpi_np=extract_mpi_np,
            extract_script=extract_script,
            extract_config_dir=extract_config_dir,
        )
        print(f"\n=== Comparison: {comparison.get('name', comparison['variable'])} ===")
        run_comparison(
            obs_dir=_path(comparison["obs_dir"]),
            model_files=model_files,
            model_labels=model_labels,
            bp_file=_path(comparison["bp_file"]),
            model_start=comparison["model_start"],
            model_time_unit=comparison.get("model_time_unit", "days"),
            missing_value=float(comparison.get("missing_value", -9999.0)),
            variable=comparison["variable"],
            ylabel=comparison.get("ylabel"),
            start=comparison.get("start"),
            end=comparison.get("end"),
            dpi=int(comparison.get("dpi", 150)),
            overview_layout=tuple(comparison.get("overview_layout", [6, 3])),
            overview_ylim=comparison.get("overview_ylim"),
            output_dir=_path(comparison["output_dir"]),
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot SCHISM station extractions against downloaded USGS data."
    )
    parser.add_argument("--config", type=Path, help="JSON comparison config")
    parser.add_argument("--comparison", action="append", help="comparison name from --config; repeatable")
    parser.add_argument("--no-extract-missing", action="store_true")
    parser.add_argument("--extract-mpi-np", type=int, default=4)
    parser.add_argument("--extract-script")
    parser.add_argument("--extract-config-dir", type=Path)
    parser.add_argument("--obs-dir", type=Path)
    parser.add_argument(
        "--model-file",
        type=Path,
        action="append",
        help="model extraction file; repeat up to three times",
    )
    parser.add_argument(
        "--model-label",
        action="append",
        help="legend label corresponding to each --model-file",
    )
    parser.add_argument("--bp-file", type=Path)
    parser.add_argument("--model-start")
    parser.add_argument("--variable")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--model-time-unit",
        choices=sorted(TIME_UNIT_TO_SECONDS),
        default="days",
    )
    parser.add_argument("--missing-value", type=float, default=-9999.0)
    parser.add_argument("--ylabel")
    parser.add_argument("--start")
    parser.add_argument("--end")
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument(
        "--overview-layout",
        type=int,
        nargs=2,
        metavar=("ROWS", "COLUMNS"),
        default=(6, 3),
        help="overview subplot layout (default: 6 3)",
    )
    parser.add_argument(
        "--overview-ylim",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="fixed y-axis range applied only to overview plots",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.config is not None:
        run_config(
            args.config,
            comparison_names=args.comparison,
            extract_missing=not args.no_extract_missing,
            extract_mpi_np=args.extract_mpi_np,
            extract_script=args.extract_script,
            extract_config_dir=args.extract_config_dir,
        )
        return

    required = {
        "--obs-dir": args.obs_dir,
        "--model-file": args.model_file,
        "--bp-file": args.bp_file,
        "--model-start": args.model_start,
        "--variable": args.variable,
        "--output-dir": args.output_dir,
    }
    missing = [name for name, value in required.items() if not value]
    if missing:
        raise ValueError(f"Missing required arguments without --config: {missing}")
    run_comparison(
        obs_dir=args.obs_dir,
        model_files=args.model_file,
        model_labels=args.model_label,
        bp_file=args.bp_file,
        model_start=args.model_start,
        model_time_unit=args.model_time_unit,
        missing_value=args.missing_value,
        variable=args.variable,
        ylabel=args.ylabel,
        start=args.start,
        end=args.end,
        dpi=args.dpi,
        overview_layout=tuple(args.overview_layout),
        overview_ylim=tuple(args.overview_ylim) if args.overview_ylim is not None else None,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
