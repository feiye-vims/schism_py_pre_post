"""Plot SCHISM elevation against CO-OPS stations for one or more runs.

The preferred JSON layout is one event file with a ``runs`` list. The first
run is treated as the main run and is plotted separately as ``ts_*.png``.
When more than one run is configured, all runs are overlaid in ``comp_*.png``.

Example::

    {
        "events": {
            "2018_v7": {
                "plot_start_day_str": "2018-09-07 00:00:00",
                "plot_end_day_str": "2018-09-22 00:00:00",
                "box": {"W": -92, "E": -88, "S": 29, "N": 31},
                "station_bp_file": "/path/to/station.in",
                "datum_shift_file": "/path/to/navd2xgeoid_shift.txt",
                "default_datum": "NAVD",
                "runs": [
                    {
                        "name": "R15a",
                        "model_start_day_str": "2018-09-07 00:00:00",
                        "elev_out_file": "/path/to/main/staout_1",
                        "shift": 0.0,
                        "line_style": "k"
                    },
                    {
                        "name": "R23",
                        "model_start_day_str": "2018-09-07 00:00:00",
                        "elev_out_file": "/path/to/other/staout_1",
                        "output_dir": "/path/to/schism/run",
                        "line_style": "--c",
                        "extract": {
                            "var_name": "elevation",
                            "bpfile": "/path/to/station.bp",
                            "nproc": 4
                        }
                    }
                ]
            }
        }
    }

The ``extract`` block is only used when ``elev_out_file`` is missing or empty.
When extracting, ``output_dir`` must contain an ``outputs`` directory with
consecutive stacks named either ``schout_1.nc``, ``schout_2.nc``, ... or
``out2d_1.nc``, ``out2d_2.nc``, ...

The legacy layout used by ``plot_coastal_act.py`` is also accepted. Additional
legacy per-run JSON files can be supplied through ``run_json_files``.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Download.download_coops_elev import get_coops_elev
from Extract.extract_schism_output import DEFAULT_EXTRACT_SCRIPT, ensure_model_file
from Grid.Bpfile import Bpfile
from Plot.plot_elev import datum_shift, get_hindcast_elev, plot_elev, write_stat


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SYMBOL_FILE = SCRIPT_DIR / "coastal_act_stats_plot_symbols.json"
DEFAULT_NEGLECT_STATIONS = ["8632837"]
DEFAULT_LINE_STYLES = ["k", "--g", "--c", "--m", "--b", "--r", ":g", ":c"]
SAMPLE_STATIONS = {"Sample": ["8720030", "8423898", "8419317", "8410140"]}
SAMPLE_DATUM = {"Sample": "NAVD"}


def ecgc_stations(station_bp_file: str) -> tuple[dict[str, list[str]], dict[str, str]]:
    noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id
    stations_groups = {
        "Florida": noaa_stations_all[:10],
        "Atlantic": noaa_stations_all[10:29],
        "GoME": noaa_stations_all[29:39],
        "GoMX_west": noaa_stations_all[41:60],
        "GoMX_east": noaa_stations_all[60:80],
        "Atlantic_inland1": noaa_stations_all[80:100],
        "Atlantic_inland2": noaa_stations_all[100:120],
        "GoMX_inland": noaa_stations_all[120:150],
        "Puerto_Rico": noaa_stations_all[150:164] + noaa_stations_all[39:41],
    }
    default_datums = {
        "Florida": "NAVD",
        "Atlantic": "NAVD",
        "GoME": "NAVD",
        "GoMX_west": "NAVD",
        "GoMX_east": "NAVD",
        "Atlantic_inland1": "NAVD",
        "Atlantic_inland2": "NAVD",
        "GoMX_inland": "NAVD",
        "Puerto_Rico": "MSL",
    }
    return stations_groups, default_datums


def subset_stations_in_box(
    box: dict[str, float],
    station_bp_file: str,
    group_name: str = "Landfall_region",
    default_datum: str = "NAVD",
    max_subplots: int = 20,
) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    station_bp = Bpfile(station_bp_file, cols=5)
    stations_all = station_bp.st_id
    stations_lonlat = station_bp.nodes[:, :2]
    idx = (
        (box["W"] - 1 < stations_lonlat[:, 0])
        * (box["E"] + 1 > stations_lonlat[:, 0])
        * (box["S"] - 1 < stations_lonlat[:, 1])
        * (box["N"] + 1 > stations_lonlat[:, 1])
    )

    stations_inbox = np.array(stations_all)[idx]
    stations_groups = {}
    default_datums = {}
    for i_group, i in enumerate(range(0, len(stations_inbox), max_subplots)):
        these_stations = stations_inbox[i:i + max_subplots]
        this_group_name = f"{group_name}_{i_group}"
        stations_groups[this_group_name] = these_stations
        default_datums[this_group_name] = default_datum
        print(f"{len(these_stations)} stations in {this_group_name}")

    return stations_groups, default_datums


def stats_scatter(
    stats: pd.DataFrame,
    var_str: str,
    region: str,
    box: dict[str, float],
    plot_symbol_dict: dict[str, Any],
    filename: str,
) -> None:
    from mpl_toolkits.basemap import Basemap

    var_str0 = var_str.split("_")[0]
    ilabel = plot_symbol_dict[region]["ilabel"]
    grid_spacing = plot_symbol_dict[region]["grid_spacing"]
    symbol_size = plot_symbol_dict[region]["symbol_size"]
    var_colorbar_str = plot_symbol_dict[var_str]["var_colorbar_str"]
    cmap = plot_symbol_dict[var_str]["cmap"]
    vmin = plot_symbol_dict[var_str]["vmin"]
    vmax = plot_symbol_dict[var_str]["vmax"]

    plt.figure(figsize=(12, 10), dpi=300)
    m = Basemap(
        projection="merc",
        resolution="f",
        llcrnrlat=box["S"] - 1,
        llcrnrlon=box["W"] - 1,
        urcrnrlat=box["N"] + 1,
        urcrnrlon=box["E"] + 1,
    )
    m.shadedrelief()
    m.drawcoastlines()
    parallels = np.arange(np.floor(box["S"]) - 1, np.ceil(box["N"]) + 1, grid_spacing)
    m.drawparallels(parallels, labels=[False, True, True, False])
    meridians = np.arange(np.ceil(box["W"]) - 1, np.ceil(box["E"]) + 1, grid_spacing)
    m.drawmeridians(meridians, labels=[False, True, True, False])
    mx, my = m(np.array(stats["station_lon"]), np.array(stats["station_lat"]))
    plt.scatter(mx, my, symbol_size, stats[var_str0].values, cmap=cmap, vmin=vmin, vmax=vmax)
    if ilabel == 1:
        for i, label in enumerate(stats["station_id"]):
            plt.annotate(f"{label}", (mx[i], my[i]))
    clb = plt.colorbar()
    clb.ax.set_title(f"{var_colorbar_str} (m)")
    plt.savefig(f"{filename}_{var_str}.png", dpi=400)
    plt.close("all")


def _load_json(json_file: str | os.PathLike[str]) -> dict[str, Any]:
    with open(json_file, "r", encoding="utf-8") as f:
        return json.load(f)


def _event_item(plot_dict: dict[str, Any], event: str) -> dict[str, Any]:
    if "events" in plot_dict:
        return plot_dict["events"][event]
    return plot_dict[event]


def _run_id(run_spec: dict[str, Any], fallback_output_dir: str | None = None) -> str:
    for key in ("name", "runid", "label"):
        if run_spec.get(key):
            return str(run_spec[key])

    output_dir = run_spec.get("output_dir", fallback_output_dir)
    if output_dir:
        runid = os.path.basename(os.path.normpath(output_dir))
        if runid:
            return runid

    elev_file = run_spec.get("elev_out_file", "model")
    return Path(elev_file).parent.name or Path(elev_file).stem


def _deduplicate_run_names(run_specs: list[dict[str, Any]]) -> None:
    seen = set()
    for i, spec in enumerate(run_specs):
        name = spec["name"]
        if name not in seen:
            seen.add(name)
            continue

        output_dir_name = ""
        if spec.get("output_dir"):
            output_dir_name = os.path.basename(os.path.normpath(spec["output_dir"]))
        candidate = output_dir_name or f"{name}_{i + 1}"
        if candidate in seen:
            candidate = f"{name}_{i + 1}"
        spec["name"] = candidate
        seen.add(candidate)


def _expand_run_specs(
    config_file: str | os.PathLike[str],
    event: str,
    run_json_files: list[str] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Return event metadata and normalized run specs.

    Event-level keys are inherited by each run. This keeps shared values such
    as ``station_bp_file`` and ``datum_shift_file`` out of every run entry.
    """

    plot_dict = _load_json(config_file)
    event_item = _event_item(plot_dict, event)
    run_items = event_item.get("runs", plot_dict.get("runs"))

    if run_items is None:
        run_items = [event_item]

    event_defaults = {k: v for k, v in event_item.items() if k != "runs"}
    if "All" in plot_dict and "box" in plot_dict["All"]:
        event_defaults["_full_domain_box"] = plot_dict["All"]["box"]
    run_specs = []
    for i, run_item in enumerate(run_items):
        spec = {**event_defaults, **run_item}
        spec["name"] = _run_id(spec, event_defaults.get("output_dir"))
        spec["shift"] = float(spec.get("shift", 0.0))
        spec["line_style"] = spec.get("line_style", DEFAULT_LINE_STYLES[i % len(DEFAULT_LINE_STYLES)])
        run_specs.append(spec)

    for json_file in run_json_files or []:
        other_dict = _load_json(json_file)
        other_item = _event_item(other_dict, event)
        spec = {**event_defaults, **other_item}
        spec["name"] = _run_id(spec, other_item.get("output_dir"))
        spec["shift"] = float(spec.get("shift", 0.0))
        spec["line_style"] = spec.get(
            "line_style", DEFAULT_LINE_STYLES[len(run_specs) % len(DEFAULT_LINE_STYLES)]
        )
        run_specs.append(spec)

    _deduplicate_run_names(run_specs)
    return event_defaults, run_specs


def _read_station_subset(run_spec: dict[str, Any]) -> range | list[int]:
    if "station_subset" in run_spec:
        return run_spec["station_subset"]

    with open(run_spec["station_bp_file"], "r", encoding="utf-8") as f:
        f.readline()
        n_station = int(f.readline().split()[0])
    return range(n_station)


def _load_model(
    run_spec: dict[str, Any],
    extract_missing: bool,
    extract_mpi_np: int,
    extract_script: str | os.PathLike[str],
    extract_config_dir: str | os.PathLike[str],
) -> pd.DataFrame:
    ensure_model_file(
        run_spec,
        extract_missing=extract_missing,
        extract_mpi_np=extract_mpi_np,
        extract_script=extract_script,
        extract_config_dir=extract_config_dir,
    )
    model = get_hindcast_elev(
        model_start_day_str=run_spec["model_start_day_str"],
        noaa_stations=None,
        station_in_file=run_spec["station_bp_file"],
        elev_out_file=run_spec["elev_out_file"],
        station_in_subset=_read_station_subset(run_spec),
    )
    model = model + run_spec.get("shift", 0.0)

    datum_shift_file = run_spec.get("datum_shift_file")
    if datum_shift_file:
        model = datum_shift(model, datum_shift_file=datum_shift_file)

    return model


def _station_groups(
    event_item: dict[str, Any],
    region: str,
    default_datum: str,
    station_bp_file: str,
    max_subplots: int,
) -> tuple[dict[str, list[str]], dict[str, str], dict[str, float]]:
    if region == "Full_domain":
        box = event_item.get("_full_domain_box") or event_item.get("box") or {"W": -100, "E": -60, "S": 8, "N": 48}
        return (*ecgc_stations(station_bp_file), box)

    if region == "FromDict":
        return (
            event_item["stations_group"],
            event_item["default_datums"],
            event_item.get("box", {"W": -100, "E": -60, "S": 8, "N": 48}),
        )

    if region == "Sample":
        return SAMPLE_STATIONS, SAMPLE_DATUM, event_item.get("box", {"W": -100, "E": -60, "S": 8, "N": 48})

    box = event_item["box"]
    stations, datums = subset_stations_in_box(
        box, station_bp_file, group_name=region, default_datum=default_datum, max_subplots=max_subplots
    )
    return stations, datums, box


def _mean_stats_line(stats: pd.DataFrame, fname: str | os.PathLike[str]) -> tuple[str, str]:
    stats_with_mean = write_stat(stats, fname)
    header = stats_with_mean.iloc[:1, :].to_string(index=False).split("\n")[0]
    values = stats_with_mean.iloc[:1, :].to_string(index=False).split("\n")[1]
    return header, values


def _plot_group_comparison(
    obs: list[pd.DataFrame],
    models: dict[str, pd.DataFrame],
    line_styles: dict[str, str],
    plot_start_day_str: str,
    plot_end_day_str: str,
    stations: list[str],
    datums: list[str],
    st_info: list[Any],
    plot_name: str,
    low_pass_filter: bool,
    nday_moving_average: int,
    subplots_shape: tuple[int | None, int | None],
) -> dict[str, pd.DataFrame]:
    stats_by_run = {}
    fig_ax = None

    for i, (runid, model) in enumerate(models.items()):
        if i == 0:
            stat, fig_ax = plot_elev(
                obs, model, plot_start_day_str, plot_end_day_str,
                stations, datums, st_info, plot_name=None, iplot=False,
                subplots_shape=subplots_shape, low_pass_filter=low_pass_filter,
                nday_moving_average=nday_moving_average,
                line_styles=["r.", line_styles[runid]],
                label_strs=["obs", runid],
            )
        else:
            stat, fig_ax = plot_elev(
                obs, model, plot_start_day_str, plot_end_day_str,
                stations, datums, st_info, plot_name=None, iplot=False,
                subplots_shape=subplots_shape, fig_ax=fig_ax,
                line_styles=[None, line_styles[runid]],
                low_pass_filter=low_pass_filter,
                nday_moving_average=nday_moving_average,
                label_strs=["obs", runid],
            )
        stats_by_run[runid] = stat

    fig_ax[0].savefig(plot_name)
    plt.close(fig_ax[0])
    return stats_by_run


def plot_coops_elev(
    config_file: str | os.PathLike[str],
    events: list[str],
    region: str = "Full_domain",
    datum: str = "NAVD",
    run_json_files: list[str] | None = None,
    output_dir: str | os.PathLike[str] = ".",
    low_pass_filter: bool = False,
    nday_moving_average: int = 0,
    outfilename_suffix: str = "",
    subplots_shape: tuple[int | None, int | None] = (10, None),
    max_subplots: int = 20,
    retrieve_method: str = "noaa_coops",
    cache_folder: str | os.PathLike[str] | None = None,
    stats_symbol_file: str | os.PathLike[str] = DEFAULT_SYMBOL_FILE,
    plot_stat: str = "MAE",
    make_scatter: bool = True,
    extract_missing: bool = True,
    extract_mpi_np: int = 4,
    extract_script: str | os.PathLike[str] = DEFAULT_EXTRACT_SCRIPT,
    extract_config_dir: str | os.PathLike[str] | None = None,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Plot CO-OPS elevation comparisons and return overall stats by event/run."""

    if not events:
        raise ValueError("events must contain at least one event name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if cache_folder is None:
        cache_folder = os.path.realpath(os.path.expanduser("~/schism10/Cache/"))
    if extract_config_dir is None:
        extract_config_dir = output_dir / "extract_configs"

    with open(stats_symbol_file, "r", encoding="utf-8") as f:
        plot_symbol_dict = json.load(f)

    all_results = {}
    for event in events:
        event_item, run_specs = _expand_run_specs(config_file, event, run_json_files=run_json_files)
        main_run = run_specs[0]
        station_bp_file = main_run["station_bp_file"]
        plot_start_day_str = event_item["plot_start_day_str"]
        plot_end_day_str = event_item["plot_end_day_str"]
        stations_groups, default_datums, box = _station_groups(
            event_item, region, datum, station_bp_file, max_subplots=max_subplots
        )
        neglect_stations = DEFAULT_NEGLECT_STATIONS + event_item.get("neglect_stations", [])

        print(f"Processing {event}: {len(run_specs)} run(s); main run is {main_run['name']}")
        models = {
            run["name"]: _load_model(
                run,
                extract_missing=extract_missing,
                extract_mpi_np=extract_mpi_np,
                extract_script=extract_script,
                extract_config_dir=extract_config_dir,
            )
            for run in run_specs
        }
        line_styles = {run["name"]: run["line_style"] for run in run_specs}
        total_stats = {run["name"]: pd.DataFrame() for run in run_specs}
        mean_stats_text = {run["name"]: "" for run in run_specs}
        final_datums = []

        for group_index, (group_name, stations) in enumerate(stations_groups.items()):
            stations = list(stations)
            group_datum = default_datums[group_name]
            filename_base = f"{event}_{group_name}_{group_datum}"
            print(f"Processing {group_name} for {event}")

            obs, datums, st_info = get_coops_elev(
                begin_time=plot_start_day_str,
                end_time=plot_end_day_str,
                noaa_stations=stations,
                default_datum=group_datum,
                cache_folder=cache_folder,
                retrieve_method=retrieve_method,
            )
            final_datums += datums

            main_stat, _ = plot_elev(
                obs, models[main_run["name"]], plot_start_day_str, plot_end_day_str,
                stations, datums, st_info,
                plot_name=str(output_dir / f"ts_{filename_base}"),
                iplot=False, subplots_shape=subplots_shape,
                line_styles=["r.", main_run["line_style"]],
                low_pass_filter=low_pass_filter,
                nday_moving_average=nday_moving_average,
                label_strs=["obs", main_run["name"]],
                figure_type="png",
            )
            total_stats[main_run["name"]] = pd.concat(
                [total_stats[main_run["name"]], main_stat], ignore_index=True
            )

            header, values = _mean_stats_line(
                main_stat, output_dir / f"stats_{main_run['name']}_{group_name}.txt"
            )
            if group_index == 0:
                mean_stats_text[main_run["name"]] += "".ljust(27) + header + "\n"
            mean_stats_text[main_run["name"]] += f"{group_name.ljust(25)}: {values}\n"

            if len(run_specs) > 1:
                comp_stats = _plot_group_comparison(
                    obs, models, line_styles, plot_start_day_str, plot_end_day_str,
                    stations, datums, st_info,
                    str(output_dir / f"comp_{filename_base}.png"),
                    low_pass_filter, nday_moving_average, subplots_shape,
                )
                for run_spec in run_specs[1:]:
                    runid = run_spec["name"]
                    total_stats[runid] = pd.concat([total_stats[runid], comp_stats[runid]], ignore_index=True)
                    header, values = _mean_stats_line(
                        comp_stats[runid], output_dir / f"stats_{runid}_{group_name}.txt"
                    )
                    if group_index == 0:
                        mean_stats_text[runid] += "".ljust(27) + header + "\n"
                    mean_stats_text[runid] += f"{group_name.ljust(25)}: {values}\n"

        suffix = f"_{outfilename_suffix}" if outfilename_suffix else ""
        overall_name = f"{event}_{region}{suffix}"
        mask = ~total_stats[main_run["name"]]["station_id"].isin(neglect_stations)

        for run_spec in run_specs:
            runid = run_spec["name"]
            run_stats = total_stats[runid]
            if run_stats.empty:
                continue
            run_mask = ~run_stats["station_id"].isin(neglect_stations)
            header, values = _mean_stats_line(
                run_stats[run_mask], output_dir / f"stats_{runid}_{overall_name}.txt"
            )
            if not mean_stats_text[runid]:
                mean_stats_text[runid] += "".ljust(27) + header + "\n"
            mean_stats_text[runid] += f"{'Overall'.ljust(25)}: {values}"
            with open(output_dir / f"mean_stats_{runid}_{overall_name}.txt", "w", encoding="utf-8") as f:
                f.write(mean_stats_text[runid])

        if make_scatter and not total_stats[main_run["name"]].empty:
            stats_scatter(
                stats=total_stats[main_run["name"]][mask],
                var_str=plot_stat,
                region=region,
                box=box,
                plot_symbol_dict=plot_symbol_dict,
                filename=str(output_dir / overall_name),
            )

        with open(output_dir / f"datum_info_{overall_name}.txt", "w", encoding="utf-8") as f:
            f.write(f"Stations with NAVD datum: {sum(np.array(final_datums) == 'NAVD')}\n")
            f.write(f"Stations with MSL datum: {sum(np.array(final_datums) == 'MSL')}\n")
            f.write(f"Stations without data: {sum(np.array(final_datums) == None)}\n")

        all_results[event] = total_stats

    return all_results


def _parse_subplot_shape(value: str) -> tuple[int | None, int | None]:
    parts = value.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("subplot shape must be ROWS,COLS, with 'none' allowed")

    parsed = []
    for part in parts:
        part = part.strip().lower()
        parsed.append(None if part in ("none", "null", "") else int(part))
    if all(x is None for x in parsed):
        raise argparse.ArgumentTypeError("at least one subplot dimension must be set")
    return tuple(parsed)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file", help="JSON config file")
    parser.add_argument("events", nargs="+", help="Event names to process")
    parser.add_argument("--run-json", action="append", default=None, help="Extra legacy per-run JSON file")
    parser.add_argument("--region", default="Full_domain")
    parser.add_argument("--datum", default="NAVD")
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--low-pass-filter", action="store_true")
    parser.add_argument("--nday-moving-average", type=int, default=0)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--subplots-shape", type=_parse_subplot_shape, default=(10, None))
    parser.add_argument("--retrieve-method", default="noaa_coops")
    parser.add_argument("--no-scatter", action="store_true")
    parser.add_argument(
        "--no-extract-missing",
        action="store_true",
        help="fail instead of extracting when a model file is missing",
    )
    parser.add_argument("--extract-mpi-np", type=int, default=4, help="MPI ranks for missing model extraction")
    parser.add_argument("--extract-script", default=str(DEFAULT_EXTRACT_SCRIPT), help="MPI extraction script")
    parser.add_argument("--extract-config-dir", default=None, help="directory for generated extraction configs")
    args = parser.parse_args()

    plot_coops_elev(
        config_file=args.config_file,
        events=args.events,
        region=args.region,
        datum=args.datum,
        run_json_files=args.run_json,
        output_dir=args.output_dir,
        low_pass_filter=args.low_pass_filter,
        nday_moving_average=args.nday_moving_average,
        outfilename_suffix=args.suffix,
        subplots_shape=args.subplots_shape,
        retrieve_method=args.retrieve_method,
        make_scatter=not args.no_scatter,
        extract_missing=not args.no_extract_missing,
        extract_mpi_np=args.extract_mpi_np,
        extract_script=args.extract_script,
        extract_config_dir=args.extract_config_dir,
    )


if __name__ == "__main__":
    main()
