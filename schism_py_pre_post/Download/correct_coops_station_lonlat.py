"""
Correct CO-OPS station lon/lat columns in a datum-shift CSV.

The intended first use is:

    python Download/correct_coops_station_lonlat.py \
        /sciclone/schism10/feiye/STOFS3D-v8/BPfiles/navd2xgeoid_shift.txt

By default the script writes a sibling file with the suffix
``.lonlat_corrected`` and an audit CSV. Use ``--in-place`` to overwrite the
input file after reviewing the audit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import time
from pathlib import Path
from typing import Any

import requests


DEFAULT_INPUT = Path(
    "/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/navd2xgeoid_shift.txt"
)
COOPS_STATION_URL = (
    "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{station}.json"
)


def _coerce_station_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Return the station dictionary from common NOAA metadata payload shapes."""
    stations = payload.get("stations")
    if isinstance(stations, list) and stations:
        return stations[0]
    if isinstance(stations, dict):
        return stations
    if "id" in payload:
        return payload
    raise ValueError(f"Unexpected NOAA metadata response shape: {payload.keys()}")


def _read_lon_lat(station: dict[str, Any]) -> tuple[float, float]:
    lon = station.get("lng", station.get("lon", station.get("longitude")))
    lat = station.get("lat", station.get("latitude"))
    if lon is None or lat is None:
        raise ValueError(f"NOAA metadata response lacks lon/lat: {station.keys()}")
    return float(lon), float(lat)


def query_coops_station_info(
    station_id: str,
    timeout: float = 20.0,
    max_attempts: int = 3,
) -> dict[str, Any]:
    """
    Query NOAA CO-OPS station metadata and normalize the core station fields.

    Returns a dict with stable fields: source, id, name, lon, lat, raw.
    """
    url = COOPS_STATION_URL.format(station=station_id)
    last_error = None
    for attempt in range(1, max_attempts + 1):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            raw = _coerce_station_payload(response.json())
            lon, lat = _read_lon_lat(raw)
            return {
                "source": "coops",
                "id": str(raw.get("id", station_id)),
                "name": raw.get("name", ""),
                "lon": lon,
                "lat": lat,
                "raw": raw,
            }
        except (requests.RequestException, json.JSONDecodeError, ValueError) as exc:
            last_error = exc
            if attempt < max_attempts:
                time.sleep(2 ** (attempt - 1))
    raise RuntimeError(f"Failed to query CO-OPS station {station_id}: {last_error}")


def approx_distance_m(lon0: float, lat0: float, lon1: float, lat1: float) -> float:
    """Small-distance approximation used only for the audit report."""
    lat_mean = math.radians((lat0 + lat1) / 2.0)
    dx = (lon1 - lon0) * 111_320.0 * math.cos(lat_mean)
    dy = (lat1 - lat0) * 110_540.0
    return math.hypot(dx, dy)


def read_shift_rows(input_file: Path) -> tuple[list[str], list[dict[str, str]]]:
    with input_file.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {input_file}")
        rows = list(reader)
    required = {"ID", "lon", "lat", "datum", "shift"}
    missing = required - set(reader.fieldnames)
    if missing:
        raise ValueError(f"{input_file} is missing required columns: {sorted(missing)}")
    return reader.fieldnames, rows


def write_shift_rows(output_file: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_audit(audit_file: Path, audit_rows: list[dict[str, Any]]) -> None:
    audit_file.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "ID",
        "datum",
        "old_lon",
        "old_lat",
        "new_lon",
        "new_lat",
        "delta_lon",
        "delta_lat",
        "approx_delta_m",
        "status",
        "note",
    ]
    with audit_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(audit_rows)


def correct_lonlat(
    input_file: Path,
    output_file: Path,
    audit_file: Path,
    precision: int = 10,
) -> tuple[int, int, int]:
    fieldnames, rows = read_shift_rows(input_file)
    station_ids = list(dict.fromkeys(row["ID"].strip() for row in rows))
    station_info: dict[str, dict[str, Any] | None] = {}

    for i, station_id in enumerate(station_ids, start=1):
        print(f"Querying CO-OPS station {i} of {len(station_ids)}: {station_id}")
        try:
            station_info[station_id] = query_coops_station_info(station_id)
        except RuntimeError as exc:
            print(f"  {exc}")
            station_info[station_id] = None

    audit_rows: list[dict[str, Any]] = []
    corrected = 0
    missing = 0

    for row in rows:
        station_id = row["ID"].strip()
        old_lon = float(row["lon"])
        old_lat = float(row["lat"])
        info = station_info[station_id]

        if info is None:
            missing += 1
            audit_rows.append({
                "ID": station_id,
                "datum": row["datum"],
                "old_lon": old_lon,
                "old_lat": old_lat,
                "new_lon": "",
                "new_lat": "",
                "delta_lon": "",
                "delta_lat": "",
                "approx_delta_m": "",
                "status": "missing_info",
                "note": "kept original lon/lat",
            })
            continue

        new_lon = float(info["lon"])
        new_lat = float(info["lat"])
        distance = approx_distance_m(old_lon, old_lat, new_lon, new_lat)
        changed = not (math.isclose(old_lon, new_lon) and math.isclose(old_lat, new_lat))

        if changed:
            corrected += 1
            row["lon"] = f"{new_lon:.{precision}f}".rstrip("0").rstrip(".")
            row["lat"] = f"{new_lat:.{precision}f}".rstrip("0").rstrip(".")

        audit_rows.append({
            "ID": station_id,
            "datum": row["datum"],
            "old_lon": old_lon,
            "old_lat": old_lat,
            "new_lon": new_lon,
            "new_lat": new_lat,
            "delta_lon": new_lon - old_lon,
            "delta_lat": new_lat - old_lat,
            "approx_delta_m": distance,
            "status": "corrected" if changed else "unchanged",
            "note": info.get("name", ""),
        })

    write_shift_rows(output_file, fieldnames, rows)
    write_audit(audit_file, audit_rows)
    return len(rows), corrected, missing


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Correct CO-OPS station lon/lat columns using NOAA metadata."
    )
    parser.add_argument("input_file", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("-o", "--output", type=Path, default=None)
    parser.add_argument("--audit", type=Path, default=None)
    parser.add_argument("--in-place", action="store_true")
    parser.add_argument("--no-backup", action="store_true")
    parser.add_argument("--precision", type=int, default=10)
    return parser.parse_args()


def next_backup_path(input_file: Path) -> Path:
    backup_file = input_file.with_suffix(input_file.suffix + ".bak_lonlat")
    if not backup_file.exists():
        return backup_file
    for i in range(1, 1000):
        candidate = input_file.with_suffix(input_file.suffix + f".bak_lonlat.{i}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Could not find an available backup name for {input_file}")


def main() -> None:
    args = parse_args()
    input_file = args.input_file

    if args.in_place:
        output_file = input_file
        write_file = input_file.with_suffix(input_file.suffix + ".tmp_lonlat_corrected")
    elif args.output is not None:
        output_file = args.output
        write_file = output_file
    else:
        output_file = input_file.with_suffix(input_file.suffix + ".lonlat_corrected")
        write_file = output_file

    audit_file = args.audit
    if audit_file is None:
        audit_file = output_file.with_suffix(output_file.suffix + ".audit.csv")

    n_rows, corrected, missing = correct_lonlat(
        input_file=input_file,
        output_file=write_file,
        audit_file=audit_file,
        precision=args.precision,
    )
    if args.in_place:
        backup_file = None
        if not args.no_backup:
            backup_file = next_backup_path(input_file)
            shutil.copy2(input_file, backup_file)
        write_file.replace(input_file)
        if backup_file is not None:
            print(f"Wrote backup file: {backup_file}")

    print(
        f"Processed {n_rows} rows; corrected {corrected}; "
        f"missing station metadata for {missing}."
    )
    print(f"Wrote corrected file: {output_file}")
    print(f"Wrote audit file: {audit_file}")


if __name__ == "__main__":
    main()
