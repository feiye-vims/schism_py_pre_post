"""Generic USGS/SCHISM time-series comparison plots.

The workflow supports any scalar variable described by :class:`VariableSpec`.
Built-in specifications are provided for salinity, temperature, and water
level.  SCHISM station extractions are expected in ``TimeHistory`` format:
the first column is model time and subsequent columns follow the station bp
file order.

USGS salinity and temperature records generally do not identify a sensor
depth.  The built-in specifications therefore label observations as
``depth not reported`` and the model as ``near-surface``.  This is an honest
screening comparison, not a depth-matched validation in stratified water.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from schism_py_pre_post.Grid.Bpfile import Bpfile
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory


@dataclass(frozen=True)
class VariableSpec:
    """Metadata and unit conversions for one model/observation variable."""

    name: str
    usgs_parameter: str
    ylabel: str
    observation_scale: float = 1.0
    model_scale: float = 1.0
    demean: bool = False
    observation_depth_label: str = "sensor depth not reported"
    model_depth_label: str = "near-surface"


VARIABLE_SPECS = {
    "salinity": VariableSpec(
        name="Salinity",
        usgs_parameter="00480",
        ylabel="Salinity (ppt/PSU)",
    ),
    "temperature": VariableSpec(
        name="Water temperature",
        usgs_parameter="00010",
        ylabel="Temperature (°C)",
    ),
    # USGS 00065 is gage height, normally reported in feet.  A datum must be
    # established separately before interpreting it as absolute elevation.
    "elevation": VariableSpec(
        name="Water level",
        usgs_parameter="00065",
        ylabel="Water level (m)",
        observation_scale=0.3048,
        observation_depth_label="gage height",
        model_depth_label="free surface",
    ),
}
VARIABLE_SPECS["gauge height"] = VARIABLE_SPECS["elevation"]


def get_variable_spec(variable: str | VariableSpec, *, demean: bool | None = None) -> VariableSpec:
    """Resolve a built-in variable name or return a user-provided spec."""

    if isinstance(variable, VariableSpec):
        spec = variable
    else:
        key = variable.strip().lower()
        if key not in VARIABLE_SPECS:
            choices = ", ".join(sorted(VARIABLE_SPECS))
            raise ValueError(f"Unknown variable {variable!r}; choose {choices} or pass VariableSpec")
        spec = VARIABLE_SPECS[key]
    return replace(spec, demean=demean) if demean is not None else spec


def _utc_naive(values) -> pd.DatetimeIndex:
    """Normalize timestamps to UTC and remove timezone metadata for plotting."""

    return pd.DatetimeIndex(pd.to_datetime(values, utc=True)).tz_localize(None)


def load_schism_station_runs(
    model_files: Mapping[str, str | Path],
    station_bp_file: str | Path,
    model_start: str,
    *,
    sec_per_time_unit: float = 1.0,
    model_loader: Callable[[str | Path], pd.DataFrame] | None = None,
) -> dict[str, pd.DataFrame]:
    """Load model runs and enforce a common station ordering.

    ``model_loader`` permits other model formats.  It must return a DataFrame
    indexed by time with station IDs as columns.
    """

    station_ids = list(Bpfile(str(station_bp_file), cols=5).make_dataframe().columns.astype(str))
    runs: dict[str, pd.DataFrame] = {}
    for run_name, model_file in model_files.items():
        if model_loader is None:
            history = TimeHistory(
                str(model_file), model_start, mask_val=-9999,
                sec_per_time_unit=sec_per_time_unit,
            )
            frame = history.df.set_index("datetime")
            if frame.shape[1] != len(station_ids):
                raise ValueError(
                    f"{model_file} has {frame.shape[1]} stations, but "
                    f"{station_bp_file} has {len(station_ids)}"
                )
            frame.columns = station_ids
        else:
            frame = model_loader(model_file).copy()
            frame.columns = frame.columns.astype(str)
            missing = [station for station in station_ids if station not in frame.columns]
            if missing:
                raise ValueError(f"Model run {run_name!r} is missing stations: {missing}")
            frame = frame.loc[:, station_ids]

        frame.index = _utc_naive(frame.index)
        frame = frame.apply(pd.to_numeric, errors="coerce").sort_index()
        runs[str(run_name)] = frame

    if not runs:
        raise ValueError("model_files must contain at least one run")
    return runs


def download_usgs_observations(
    station_ids: Sequence[str],
    spec: VariableSpec,
    start: str | pd.Timestamp,
    end: str | pd.Timestamp,
    *,
    cache_file: str | Path | None = None,
) -> dict[str, object]:
    """Download USGS instantaneous values and key them by requested site ID."""

    # Import lazily so model loading and plotting can be used without the
    # downloader's optional dependencies.
    from schism_py_pre_post.Download.download_usgs import download_stations

    requested = [str(station) for station in station_ids]
    data = download_stations(
        param_id=spec.usgs_parameter,
        station_ids=requested,
        cache_fname=None if cache_file is None else str(cache_file),
        datelist=pd.date_range(start=start, end=end),
    )
    downloaded = {str(item.station_info["id"]): item for item in data}
    return {station: downloaded.get(station) for station in requested}


def _observation_series(item: object, scale: float) -> pd.Series:
    if item is None or item.df.empty:
        return pd.Series(dtype=float, index=pd.DatetimeIndex([]))
    values = pd.to_numeric(item.df["value"], errors="coerce").to_numpy(dtype=float) * scale
    series = pd.Series(values, index=_utc_naive(item.df["date"]))
    return series[~series.index.duplicated(keep="last")].sort_index().dropna()


def plot_station_comparisons(
    model_runs: Mapping[str, pd.DataFrame],
    observations: Mapping[str, object],
    spec: VariableSpec,
    start: str | pd.Timestamp,
    end: str | pd.Timestamp,
    output_dir: str | Path,
    *,
    output_prefix: str = "usgs",
    output_format: str = "png",
    stations_per_figure: int = 12,
    columns: int = 3,
    dpi: int = 180,
) -> list[Path]:
    """Plot observation/model comparisons and return the written paths."""

    if not model_runs:
        raise ValueError("model_runs must contain at least one run")
    if stations_per_figure < 1 or columns < 1:
        raise ValueError("stations_per_figure and columns must be positive")
    first_run = next(iter(model_runs.values()))
    station_ids = list(first_run.columns.astype(str))
    if not station_ids:
        raise ValueError("model runs must contain at least one station")
    for run_name, frame in model_runs.items():
        if list(frame.columns.astype(str)) != station_ids:
            raise ValueError(f"Model run {run_name!r} has a different station order")

    start_time, end_time = _utc_naive([start, end])
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = int(np.ceil(min(stations_per_figure, len(station_ids)) / columns))
    written: list[Path] = []

    for chunk_number, offset in enumerate(range(0, len(station_ids), stations_per_figure), start=1):
        chunk = station_ids[offset:offset + stations_per_figure]
        fig, axes = plt.subplots(
            rows, columns, figsize=(6.3 * columns, 3.5 * rows),
            squeeze=False, sharex=True,
        )
        flat_axes = axes.ravel()
        for axis, station_id in zip(flat_axes, chunk):
            obs = _observation_series(observations.get(station_id), spec.observation_scale)
            if not obs.empty:
                obs = obs.loc[(obs.index >= start_time) & (obs.index <= end_time)]
            if spec.demean and not obs.empty:
                obs = obs - obs.mean()
            if not obs.empty:
                axis.plot(obs.index, obs, ".", color="tab:red", markersize=3, label="USGS")

            for run_name, frame in model_runs.items():
                model = frame[station_id] * spec.model_scale
                model = model.loc[(model.index >= start_time) & (model.index <= end_time)]
                if spec.demean and not model.empty:
                    model = model - model.mean()
                axis.plot(model.index, model, linewidth=1.2, label=str(run_name))

            axis.set_title(station_id)
            axis.grid(alpha=0.3)
            axis.tick_params(axis="x", labelrotation=25)
            if axis is flat_axes[0]:
                axis.legend(fontsize="small")
            if axis.get_subplotspec().is_first_col():
                axis.set_ylabel(spec.ylabel + ("; demeaned" if spec.demean else ""))

        for axis in flat_axes[len(chunk):]:
            axis.set_visible(False)

        fig.suptitle(
            f"{spec.name}: USGS ({spec.observation_depth_label}) vs "
            f"model ({spec.model_depth_label})",
            fontsize=14,
        )
        fig.tight_layout()
        filename = (
            f"{output_prefix}_{spec.name.lower().replace(' ', '_')}_"
            f"{chunk_number:02d}.{output_format}"
        )
        path = output_dir / filename
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    return written


def plot_usgs_variable(
    variable: str | VariableSpec,
    station_bp_file: str | Path,
    model_start: str,
    plot_start: str,
    plot_end: str,
    output_dir: str | Path,
    model_files: Mapping[str, str | Path],
    *,
    sec_per_time_unit: float = 1.0,
    cache_dir: str | Path | None = None,
    output_prefix: str = "usgs",
    demean: bool | None = None,
    model_loader: Callable[[str | Path], pd.DataFrame] | None = None,
    **plot_options,
) -> list[Path]:
    """Run the complete USGS download, SCHISM load, and plot workflow.

    Example::

        plot_usgs_variable(
            "salinity", "stations.bp", "2024-03-05 00:00:00",
            "2024-03-10 00:00:00", "2024-04-10 00:00:00", "plots",
            {"run-a": "salinity.stations.dat"},
            sec_per_time_unit=86400, cache_dir="plots/cache",
        )
    """

    spec = get_variable_spec(variable, demean=demean)
    runs = load_schism_station_runs(
        model_files, station_bp_file, model_start,
        sec_per_time_unit=sec_per_time_unit, model_loader=model_loader,
    )
    station_ids = list(next(iter(runs.values())).columns.astype(str))

    cache_file = None
    if cache_dir is not None:
        safe_dates = (
            f"{pd.Timestamp(plot_start):%Y%m%d%H%M}-"
            f"{pd.Timestamp(plot_end):%Y%m%d%H%M}"
        )
        cache_file = Path(cache_dir) / (
            f"usgs_{spec.usgs_parameter}_{Path(station_bp_file).stem}_{safe_dates}.csv"
        )
        cache_file.parent.mkdir(parents=True, exist_ok=True)

    observations = download_usgs_observations(
        station_ids, spec, plot_start, plot_end, cache_file=cache_file,
    )
    return plot_station_comparisons(
        runs, observations, spec, plot_start, plot_end, output_dir,
        output_prefix=output_prefix, **plot_options,
    )
