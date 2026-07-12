"""
Download USGS continuous observations for every matching station in a polygon.

The output directory is resumable: each station/parameter request is written as
an independent Parquet or CSV part, and rerunning an identical request skips
completed parts.
GeoJSON polygons require no geospatial dependency. Other vector formats use
GeoPandas when it is installed.

Examples:

    python Download/download_usgs_polygon.py domain.geojson \
        --start 2024-01-01 --end 2024-01-02 \
        --parameter streamflow --output-dir output/usgs_streamflow

    USGS_API_KEY=... python Download/download_usgs_polygon.py domain.shp \
        --start 2024-01-01 --end 2024-02-01 \
        --parameter salinity conductance --derive salinity \
        --output-dir output/usgs_salinity

    python Download/download_usgs_polygon.py domain.shp \
        --start 2024-01-01 --end 2024-02-01 \
        --parameter salinity_best_available \
        --output-dir output/usgs_salinity_best

Common parameter aliases and output units are streamflow (m^3/s), gage_height
(m), temperature (degC), salinity (PSU), conductance (mS/cm),
dissolved_oxygen (mg/L), ph (1), and turbidity (FNU). All observation times
are timezone-aware UTC.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

try:
    from .download_obs import (
        ObsRequest,
        USGS_PARAMETER_ALIASES,
        download_obs,
        feature_lon_lat,
        iso_interval,
        iter_ogc_features,
        make_session,
        parse_time,
    )
except ImportError:
    from download_obs import (
        ObsRequest,
        USGS_PARAMETER_ALIASES,
        download_obs,
        feature_lon_lat,
        iso_interval,
        iter_ogc_features,
        make_session,
        parse_time,
    )


Polygon = list[list[tuple[float, float]]]

SALINITY_BEST_AVAILABLE = "salinity_best_available"
SALINITY_BEST_AVAILABLE_ALIASES = {
    SALINITY_BEST_AVAILABLE,
    "salinity_generic",
}


def _normalize_requested_parameter(value: str) -> str:
    normalized = value.lower()
    if normalized in SALINITY_BEST_AVAILABLE_ALIASES:
        return SALINITY_BEST_AVAILABLE
    return USGS_PARAMETER_ALIASES.get(normalized, value)


def _normalize_parameter_code(value: Any) -> str:
    parameter = str(value)
    return parameter.zfill(5) if parameter.isdigit() else parameter


def _geometry_polygons(geometry: dict[str, Any] | None) -> list[Polygon]:
    if not geometry:
        return []
    geometry_type = geometry.get("type")
    coordinates = geometry.get("coordinates", [])
    if geometry_type == "Polygon":
        coordinates = [coordinates]
    elif geometry_type != "MultiPolygon":
        raise ValueError(f"Expected Polygon or MultiPolygon, got {geometry_type!r}")
    return [
        [[(float(x), float(y)) for x, y, *_ in ring] for ring in polygon]
        for polygon in coordinates
    ]


def _geojson_polygons(payload: dict[str, Any]) -> list[Polygon]:
    payload_type = payload.get("type")
    if payload_type == "FeatureCollection":
        polygons = []
        for feature in payload.get("features", []):
            polygons.extend(_geometry_polygons(feature.get("geometry")))
        return polygons
    if payload_type == "Feature":
        return _geometry_polygons(payload.get("geometry"))
    return _geometry_polygons(payload)


def read_polygons(path: Path) -> list[Polygon]:
    if path.suffix.lower() in {".json", ".geojson"}:
        with path.open(encoding="utf-8") as stream:
            polygons = _geojson_polygons(json.load(stream))
    else:
        try:
            import geopandas as gpd
        except ImportError as exc:
            raise RuntimeError(
                "GeoPandas is required for polygon formats other than GeoJSON"
            ) from exc
        frame = gpd.read_file(path)
        if frame.crs is None:
            raise ValueError(f"Polygon file has no CRS: {path}")
        frame = frame.to_crs("EPSG:4326")
        polygons = []
        for geometry in frame.geometry:
            if geometry is not None and not geometry.is_empty:
                polygons.extend(_geometry_polygons(geometry.__geo_interface__))
    if not polygons:
        raise ValueError(f"No Polygon or MultiPolygon geometry found in {path}")
    return polygons


def polygon_bounds(polygons: list[Polygon]) -> tuple[float, float, float, float]:
    points = [point for polygon in polygons for ring in polygon for point in ring]
    x, y = zip(*points)
    return min(x), min(y), max(x), max(y)


def _point_on_segment(
    point: tuple[float, float],
    start: tuple[float, float],
    end: tuple[float, float],
    tolerance: float = 1e-12,
) -> bool:
    x, y = point
    x1, y1 = start
    x2, y2 = end
    cross = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
    if abs(cross) > tolerance:
        return False
    return (
        min(x1, x2) - tolerance <= x <= max(x1, x2) + tolerance
        and min(y1, y2) - tolerance <= y <= max(y1, y2) + tolerance
    )


def _ring_location(point: tuple[float, float], ring: list[tuple[float, float]]) -> int:
    """Return 1 inside, 0 outside, and 2 on the ring boundary."""
    inside = False
    x, y = point
    for start, end in zip(ring, ring[1:] + ring[:1]):
        if _point_on_segment(point, start, end):
            return 2
        x1, y1 = start
        x2, y2 = end
        if (y1 > y) != (y2 > y):
            intersection_x = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < intersection_x:
                inside = not inside
    return 1 if inside else 0


def polygon_covers(point: tuple[float, float], polygon: Polygon) -> bool:
    outer = _ring_location(point, polygon[0])
    if outer == 0:
        return False
    for hole in polygon[1:]:
        location = _ring_location(point, hole)
        if location == 1:
            return False
        if location == 2:
            return True
    return True


def polygons_cover(point: tuple[float, float], polygons: list[Polygon]) -> bool:
    return any(polygon_covers(point, polygon) for polygon in polygons)


def _series_overlaps(properties: dict[str, Any], start: Any, end: Any) -> bool:
    series_start = parse_time(properties.get("begin_utc") or properties.get("begin"))
    series_end = parse_time(properties.get("end_utc") or properties.get("end"))
    request_start = parse_time(start)
    request_end = parse_time(end)
    return not (
        series_start is not None and request_end is not None and series_start > request_end
        or series_end is not None and request_start is not None and series_end < request_start
    )


def discover_usgs_series(
    polygons: list[Polygon],
    parameter: str,
    start: Any,
    end: Any,
    session: Any,
) -> pd.DataFrame:
    bbox = ",".join(f"{value:.12g}" for value in polygon_bounds(polygons))
    features = iter_ogc_features(
        session,
        "time-series-metadata",
        {"bbox": bbox, "parameter_code": parameter},
    )
    rows = []
    for feature in features:
        properties = feature.get("properties", {})
        station_id = properties.get("monitoring_location_id")
        lon, lat = feature_lon_lat(feature)
        if (
            not station_id
            or lon is None
            or lat is None
            or not polygons_cover((lon, lat), polygons)
        ):
            continue
        if str(properties.get("computation_identifier", "")).lower() != "instantaneous":
            continue
        if not _series_overlaps(properties, start, end):
            continue
        rows.append({
            "station_id": station_id,
            "time_series_id": feature.get("id"),
            "lon": lon,
            "lat": lat,
            "parameter": properties.get("parameter_code") or parameter,
            "unit": properties.get("unit_of_measure"),
            "primary": properties.get("primary"),
            "begin": properties.get("begin_utc") or properties.get("begin"),
            "end": properties.get("end_utc") or properties.get("end"),
        })
    series = pd.DataFrame(rows)
    if series.empty:
        return series

    # Prefer reviewed primary series when USGS publishes both primary and
    # short-lived non-primary series for the same station and parameter.
    primary = series["primary"].fillna("").astype(str).str.lower() == "primary"
    stations_with_primary = set(series.loc[primary, "station_id"])
    series = series[primary | ~series["station_id"].isin(stations_with_primary)]
    return series.drop_duplicates("time_series_id").sort_values(
        ["station_id", "time_series_id"]
    ).reset_index(drop=True)


def select_requested_series(
    series: pd.DataFrame, parameters: list[str]
) -> pd.DataFrame:
    if series.empty:
        return series
    series = series.copy()
    series["parameter"] = series["parameter"].map(_normalize_parameter_code)
    selected = []
    for parameter in parameters:
        if parameter != SALINITY_BEST_AVAILABLE:
            selected.append(series[series["parameter"] == parameter])
            continue
        salinity = series[series["parameter"].isin(["00480", "00095"])]
        observed_stations = set(
            salinity.loc[salinity["parameter"] == "00480", "station_id"]
        )
        selected.append(salinity[
            (salinity["parameter"] == "00480")
            | ~salinity["station_id"].isin(observed_stations)
        ])
    if not selected:
        return series.iloc[0:0]
    return (
        pd.concat(selected, ignore_index=True)
        .drop_duplicates("time_series_id")
        .sort_values(["station_id", "parameter", "time_series_id"])
        .reset_index(drop=True)
    )


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(temporary, path)


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def _write_part(frame: pd.DataFrame, path: Path, output_format: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    if output_format == "parquet":
        frame.to_parquet(temporary, index=False, compression="zstd")
    else:
        frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def _read_part(path: Path, output_format: str) -> pd.DataFrame:
    return pd.read_parquet(path) if output_format == "parquet" else pd.read_csv(path)


def _safe_station_name(station_id: str) -> str:
    return station_id.replace("/", "_").replace("\\", "_")


def _request_manifest(
    args: argparse.Namespace, parameters: list[str]
) -> dict[str, Any]:
    polygon_digest = hashlib.sha256(args.polygon.read_bytes()).hexdigest()
    return {
        "version": 4,
        "polygon": str(args.polygon.resolve()),
        "polygon_sha256": polygon_digest,
        "start": iso_interval(args.start, args.end).split("/", 1)[0],
        "end": iso_interval(args.start, args.end).split("/", 1)[1],
        "parameters": parameters,
        "derive": args.derive,
        "format": args.format,
        "time_zone": "UTC",
        "standard_units": True,
    }


def _check_or_create_run(output_dir: Path, manifest: dict[str, Any]) -> None:
    manifest_path = output_dir / "manifest.json"
    output_dir.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8") as stream:
            previous = json.load(stream)
        previous_request = {key: previous.get(key) for key in manifest}
        if previous_request != manifest:
            raise RuntimeError(
                f"Output directory belongs to a different request: {output_dir}"
            )
    elif any(output_dir.iterdir()):
        raise RuntimeError(f"Refusing to use non-empty output directory: {output_dir}")
    else:
        _atomic_json(manifest, manifest_path)


def _write_mean_xyz(parts: Iterable[Path], output_format: str, path: Path) -> None:
    rows = []
    for part in parts:
        observations = _read_part(part, output_format)
        values = pd.to_numeric(observations.get("value"), errors="coerce")
        if observations.empty or not values.notna().any():
            continue
        derived_from = observations.get("derived_from_parameter")
        is_conductance_derived = (
            derived_from is not None
            and pd.to_numeric(derived_from, errors="coerce").eq(95).any()
        )
        rows.append({
            "lon": observations["lon"].iloc[0],
            "lat": observations["lat"].iloc[0],
            "mean_value": values.mean(),
            "station_id": observations["station_id"].iloc[0],
            "station_name": (
                observations["station_name"].iloc[0]
                if "station_name" in observations
                else ""
            ),
            "variable": observations["variable"].iloc[0],
            "unit": observations["unit"].iloc[0],
            "derivation": (
                "derived_from_conductance" if is_conductance_derived else "observed"
            ),
        })
    columns = [
        "lon", "lat", "mean_value", "station_id", "station_name", "variable",
        "unit", "derivation"
    ]
    mean = pd.DataFrame(rows, columns=columns)
    temporary = path.with_suffix(path.suffix + ".tmp")
    mean[["lon", "lat", "mean_value", "station_id", "derivation"]].to_csv(
        temporary, sep=" ", index=False, header=False
    )
    os.replace(temporary, path)
    _atomic_csv(mean, path.with_suffix(".csv"))

    bp_path = path.with_suffix(".bp")
    bp_temporary = bp_path.with_suffix(bp_path.suffix + ".tmp")
    with bp_temporary.open("w", encoding="utf-8") as stream:
        stream.write(f"{bp_path.name}\n")
        stream.write(f"{len(mean.index)}\n")
        for bp_id, row in enumerate(mean.itertuples(index=False), start=1):
            observation_type = (
                "derived" if row.derivation != "observed" else "observed"
            )
            stream.write(
                f"{bp_id} {row.lon:.12g} {row.lat:.12g} 0 "
                f"! station_id={json.dumps(str(row.station_id))} "
                f"station_name={json.dumps(str(row.station_name))} "
                f"variable={json.dumps(str(row.variable))} "
                f"observation_type={observation_type} "
                f"derivation={row.derivation}\n"
            )
    os.replace(bp_temporary, bp_path)


def run(args: argparse.Namespace) -> None:
    requested_parameters = (
        args.parameter if isinstance(args.parameter, list) else [args.parameter]
    )
    parameters = list(dict.fromkeys(
        _normalize_requested_parameter(value) for value in requested_parameters
    ))
    best_available_requested = SALINITY_BEST_AVAILABLE in parameters
    if best_available_requested and any(
        parameter in {"00480", "00095"} for parameter in parameters
    ):
        raise ValueError(
            "salinity_best_available cannot be combined with explicit salinity "
            "or conductance"
        )
    if args.derive and "00095" not in parameters and not best_available_requested:
        raise ValueError(
            "--derive salinity requires parameter 00095 or conductance"
        )
    discovery_parameters = list(dict.fromkeys(
        discovered_parameter
        for parameter in parameters
        for discovered_parameter in (
            ["00480", "00095"]
            if parameter == SALINITY_BEST_AVAILABLE
            else [parameter]
        )
    ))

    polygons = read_polygons(args.polygon)
    manifest = _request_manifest(args, parameters)
    _check_or_create_run(args.output_dir, manifest)
    parts_dir = args.output_dir / "observations"
    parts_dir.mkdir(exist_ok=True)

    session = make_session()
    api_key = os.environ.get("USGS_API_KEY")
    if api_key:
        session.headers["X-Api-Key"] = api_key

    discovery_path = args.output_dir / "discovered_series.csv"
    if discovery_path.exists() and not args.refresh_discovery:
        discovered_series = pd.read_csv(
            discovery_path,
            dtype={"station_id": str, "time_series_id": str, "parameter": str},
        )
        print(f"Using cached station discovery: {discovery_path}")
    else:
        discovered = [
            discover_usgs_series(
                polygons, parameter, args.start, args.end, session=session
            )
            for parameter in discovery_parameters
        ]
        discovered_series = pd.concat(discovered, ignore_index=True)
        _atomic_csv(discovered_series, discovery_path)
    series = select_requested_series(discovered_series, parameters)
    _atomic_csv(series, args.output_dir / "selected_series.csv")
    if series.empty:
        print("No matching instantaneous USGS time series found.")
        return

    selected_series = defaultdict(set)
    for row in series.itertuples(index=False):
        parameter = _normalize_parameter_code(row.parameter)
        selected_series[(row.station_id, parameter)].add(row.time_series_id)

    failures_path = args.output_dir / "failures.csv"
    failures = (
        pd.read_csv(
            failures_path, dtype={"station_id": str, "parameter": str}
        )
        if failures_path.exists()
        else pd.DataFrame(columns=["station_id", "parameter", "error"])
    )
    requests_to_download = sorted(selected_series)
    extension = ".parquet" if args.format == "parquet" else ".csv"
    stations_path = args.output_dir / "stations.csv"
    stations_table = (
        pd.read_csv(stations_path, dtype={"station_id": str})
        if stations_path.exists()
        else pd.DataFrame()
    )
    for number, (station_id, parameter) in enumerate(
        requests_to_download, start=1
    ):
        part = parts_dir / (
            f"{_safe_station_name(station_id)}__{parameter}{extension}"
        )
        failure_matches = (
            failures["station_id"].eq(station_id)
            & failures["parameter"].eq(parameter)
        )
        if part.exists():
            print(
                f"[{number}/{len(requests_to_download)}] "
                f"{station_id} {parameter}: already complete"
            )
            failures = failures[~failure_matches]
            continue
        print(
            f"[{number}/{len(requests_to_download)}] "
            f"{station_id} {parameter}: downloading"
        )
        try:
            derivation = (
                "salinity"
                if parameter == "00095" and best_available_requested
                else args.derive if parameter == "00095" else None
            )
            observations, stations = download_obs(
                ObsRequest(
                    source="usgs",
                    station_ids=[station_id],
                    start=args.start,
                    end=args.end,
                    parameter=parameter,
                    derive=derivation,
                ),
                session=session,
            )
            if "time_series_id" in observations:
                observations = observations[
                    observations["time_series_id"].isin(
                        selected_series[(station_id, parameter)]
                    )
                ].reset_index(drop=True)
            if derivation:
                observations["observation_type"] = "derived"
                observations["derivation"] = "derived_from_conductance"
            else:
                observations["observation_type"] = "observed"
                observations["derivation"] = "observed"
            _write_part(observations, part, args.format)
            stations_table = pd.concat(
                [stations_table, stations], ignore_index=True
            ).drop_duplicates("station_id", keep="last")
            _atomic_csv(stations_table, stations_path)
            failures = failures[~failure_matches]
            _atomic_csv(failures, failures_path)
        except Exception as exc:
            failure = pd.DataFrame([{
                "station_id": station_id,
                "parameter": parameter,
                "error": str(exc),
            }])
            failures = pd.concat(
                [failures[~failure_matches], failure],
                ignore_index=True,
            )
            _atomic_csv(failures, failures_path)
            if args.fail_fast:
                raise

    _atomic_csv(failures, failures_path)
    completed_parts = sorted(parts_dir.glob(f"*{extension}"))
    _write_mean_xyz(completed_parts, args.format, args.output_dir / "mean_xyz.txt")
    station_ids = {station_id for station_id, _ in requests_to_download}
    manifest.update({
        "selected_station_count": len(station_ids),
        "selected_station_parameter_count": len(requests_to_download),
        "completed_station_parameter_count": len(completed_parts),
        "failed_station_parameter_count": len(failures.index),
    })
    _atomic_json(manifest, args.output_dir / "manifest.json")
    print(
        f"Completed {len(completed_parts)} of {len(requests_to_download)} "
        f"station/parameter requests; "
        f"{len(failures.index)} failed. Output: {args.output_dir}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download USGS continuous observations within a polygon."
    )
    parser.add_argument("polygon", type=Path)
    parser.add_argument("--start", required=True)
    parser.add_argument("--end", required=True)
    parser.add_argument(
        "--parameter",
        nargs="+",
        required=True,
        help=(
            "one or more aliases or USGS parameter codes; "
            "salinity_best_available prefers observed salinity at each station "
            "and otherwise derives it from conductance"
        ),
    )
    parser.add_argument("--derive", choices=["salinity"])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--format", choices=["parquet", "csv"], default="parquet")
    parser.add_argument(
        "--refresh-discovery",
        action="store_true",
        help="query USGS again instead of reusing selected_series.csv",
    )
    parser.add_argument("--fail-fast", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
