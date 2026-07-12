"""
Generic observation downloader for CO-OPS, USGS Water Data, and NDBC.

The public return shape is intentionally source-neutral:

``download_obs(...)`` returns ``(obs, stations)`` where ``obs`` is a tidy
long-format DataFrame with columns like:

    source, station_id, station_name, lon, lat, time, variable, parameter,
    value, unit, quality, datum

and ``stations`` contains one normalized station-info row per station.

USGS uses the modern Water Data OGC API at https://api.waterdata.usgs.gov/.
The older NWIS/RDB services are deliberately not used here.

Sample CLI commands:

    python Download/download_obs.py coops 8725520 \
        --start 2024-01-01 --end 2024-01-02 \
        --variable water_level --datum MSL

    python Download/download_obs.py usgs 07374000 \
        --start 2024-01-01 --end 2024-01-01T01:00:00Z \
        --parameter streamflow

    python Download/download_obs.py usgs 07264000 \
        --start 2024-01-01 --end 2024-01-01T01:00:00Z \
        --parameter conductance --derive salinity

    python Download/download_obs.py ndbc 44025 --variable WVHT
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


COOPS_DATA_URL = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
COOPS_STATION_URL = "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{station}.json"
USGS_COLLECTION_URL = "https://api.waterdata.usgs.gov/ogcapi/v0/collections/{collection}/items"
NDBC_REALTIME_URL = "https://www.ndbc.noaa.gov/data/realtime2/{station}.txt"
NDBC_HISTORICAL_URL = "https://www.ndbc.noaa.gov/data/historical/stdmet/{station}h{year}.txt.gz"
NDBC_STATION_PAGE_URL = "https://www.ndbc.noaa.gov/station_page.php?station={station}"


USGS_PARAMETER_ALIASES = {
    "streamflow": "00060",
    "discharge": "00060",
    "gage_height": "00065",
    "gage height": "00065",
    "gauge_height": "00065",
    "gauge height": "00065",
    "temperature": "00010",
    "water_temperature": "00010",
    "water temperature": "00010",
    "conductance": "00095",
    "specific_conductance": "00095",
    "specific conductance": "00095",
    "salinity": "00480",
    "dissolved_oxygen": "00300",
    "dissolved oxygen": "00300",
    "do": "00300",
    "ph": "00400",
    "turbidity": "63680",
}

USGS_STANDARD_VARIABLES = {
    "00060": "streamflow",
    "00065": "gage_height",
    "00010": "water_temperature",
    "00095": "specific_conductance",
    "00480": "salinity",
    "00300": "dissolved_oxygen",
    "00400": "ph",
    "63680": "turbidity",
}

# Keys are normalized by _normalize_unit(). Conversion is target = value *
# scale + offset. Direct USGS salinity is reported in parts per thousand; its
# numeric value is retained as the conventional practical-salinity product.
USGS_STANDARD_UNITS = {
    "00060": {
        "ft^3/s": ("m^3/s", 0.028316846592, 0.0),
        "ft3/s": ("m^3/s", 0.028316846592, 0.0),
        "cfs": ("m^3/s", 0.028316846592, 0.0),
        "m^3/s": ("m^3/s", 1.0, 0.0),
        "m3/s": ("m^3/s", 1.0, 0.0),
    },
    "00065": {
        "ft": ("m", 0.3048, 0.0),
        "feet": ("m", 0.3048, 0.0),
        "m": ("m", 1.0, 0.0),
    },
    "00010": {
        "degc": ("degC", 1.0, 0.0),
        "deg c": ("degC", 1.0, 0.0),
        "degree_celsius": ("degC", 1.0, 0.0),
        "degrees celsius": ("degC", 1.0, 0.0),
        "degf": ("degC", 5.0 / 9.0, -32.0 * 5.0 / 9.0),
        "deg f": ("degC", 5.0 / 9.0, -32.0 * 5.0 / 9.0),
    },
    "00095": {
        "us/cm": ("mS/cm", 0.001, 0.0),
        "umho/cm": ("mS/cm", 0.001, 0.0),
        "microsiemens per centimeter": ("mS/cm", 0.001, 0.0),
        "ms/cm": ("mS/cm", 1.0, 0.0),
        "s/m": ("mS/cm", 10.0, 0.0),
    },
    "00480": {
        "ppth": ("PSU", 1.0, 0.0),
        "ppt": ("PSU", 1.0, 0.0),
        "parts per thousand": ("PSU", 1.0, 0.0),
        "psu": ("PSU", 1.0, 0.0),
        "1": ("PSU", 1.0, 0.0),
    },
    "00300": {
        "mg/l": ("mg/L", 1.0, 0.0),
    },
    "00400": {
        "std units": ("1", 1.0, 0.0),
        "std": ("1", 1.0, 0.0),
        "standard units": ("1", 1.0, 0.0),
        "1": ("1", 1.0, 0.0),
    },
    "63680": {
        "fnu": ("FNU", 1.0, 0.0),
    },
}


@dataclass(frozen=True)
class ObsRequest:
    source: str
    station_ids: list[str]
    start: Any | None = None
    end: Any | None = None
    variable: str | None = None
    parameter: str | None = None
    derive: str | None = None
    datum: str | None = None
    units: str = "metric"


def make_session() -> requests.Session:
    retry = Retry(
        total=5,
        connect=3,
        read=3,
        status=3,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET"},
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=16, pool_maxsize=16)
    session = requests.Session()
    session.mount("https://", adapter)
    session.headers.update({
        "User-Agent": "schism-py-pre-post observation downloader",
        "Accept": "application/json, text/csv, text/plain, */*",
    })
    return session


def normalize_station_id(source: str, station_id: str) -> str:
    station_id = str(station_id).strip()
    if source.lower() == "usgs" and not station_id.upper().startswith("USGS-"):
        return f"USGS-{station_id}"
    return station_id


def parse_time(value: Any | None) -> pd.Timestamp | None:
    if value is None:
        return None
    return pd.to_datetime(value, utc=True)


def iso_interval(start: Any | None, end: Any | None) -> str | None:
    start_ts = parse_time(start)
    end_ts = parse_time(end)
    if start_ts is None and end_ts is None:
        return None
    if start_ts is None or end_ts is None:
        raise ValueError("Both start and end are required for interval queries")
    start_s = start_ts.isoformat().replace("+00:00", "Z")
    end_s = end_ts.isoformat().replace("+00:00", "Z")
    return f"{start_s}/{end_s}"


def coops_date(value: Any) -> str:
    return pd.to_datetime(value).strftime("%Y%m%d")


def feature_lon_lat(feature: dict[str, Any]) -> tuple[float | None, float | None]:
    coords = feature.get("geometry", {}).get("coordinates")
    if isinstance(coords, list) and len(coords) >= 2:
        return float(coords[0]), float(coords[1])
    return None, None


def request_json(session: requests.Session, url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    response = session.get(url, params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def iter_ogc_features(
    session: requests.Session,
    collection: str,
    params: dict[str, Any],
) -> list[dict[str, Any]]:
    url = USGS_COLLECTION_URL.format(collection=collection)
    params = {"f": "json", "limit": 10000, **params}
    features: list[dict[str, Any]] = []

    while url:
        payload = request_json(session, url, params=params)
        features.extend(payload.get("features", []))
        next_url = None
        for link in payload.get("links", []):
            if link.get("rel") == "next" and link.get("href"):
                next_url = urljoin(url, link["href"])
                break
        url = next_url
        params = None

    return features


def get_usgs_station_info(
    station_ids: list[str],
    session: requests.Session | None = None,
) -> pd.DataFrame:
    session = session or make_session()
    rows = []
    for station_id in station_ids:
        mlid = normalize_station_id("usgs", station_id)
        features = iter_ogc_features(
            session,
            "monitoring-locations",
            {"id": mlid},
        )
        if not features:
            rows.append({"source": "usgs", "station_id": mlid, "status": "missing"})
            continue
        feature = features[0]
        props = feature.get("properties", {})
        lon, lat = feature_lon_lat(feature)
        rows.append({
            "source": "usgs",
            "station_id": props.get("id", mlid),
            "station_name": props.get("monitoring_location_name"),
            "lon": lon,
            "lat": lat,
            "vertical_datum": props.get("vertical_datum"),
            "status": "ok",
        })
    return pd.DataFrame(rows)


def download_usgs_continuous(
    station_ids: list[str],
    start: Any,
    end: Any,
    parameter: str = "00060",
    derive: str | None = None,
    session: requests.Session | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    session = session or make_session()
    parameter = USGS_PARAMETER_ALIASES.get(str(parameter).lower(), str(parameter))
    derive = derive.lower() if derive is not None else None
    if derive not in (None, "salinity", "salinity_from_conductance"):
        raise ValueError("USGS derive must be None or 'salinity'")
    if derive in ("salinity", "salinity_from_conductance") and parameter != "00095":
        raise ValueError("Salinity derivation requires USGS specific conductance parameter 00095")
    stations = get_usgs_station_info(station_ids, session=session)
    station_map = stations.set_index("station_id").to_dict("index") if not stations.empty else {}
    rows = []

    for station_id in station_ids:
        mlid = normalize_station_id("usgs", station_id)
        features = iter_ogc_features(
            session,
            "continuous",
            {
                "monitoring_location_id": mlid,
                "parameter_code": parameter,
                "datetime": iso_interval(start, end),
            },
        )
        info = station_map.get(mlid, {})
        for feature in features:
            props = feature.get("properties", {})
            lon, lat = feature_lon_lat(feature)
            rows.append({
                "source": "usgs",
                "station_id": props.get("monitoring_location_id", mlid),
                "station_name": props.get("monitoring_location_name") or info.get("station_name"),
                "lon": lon if lon is not None else info.get("lon"),
                "lat": lat if lat is not None else info.get("lat"),
                "time": props.get("time"),
                "variable": props.get("parameter_code", parameter),
                "parameter": props.get("parameter_code", parameter),
                "value": props.get("value"),
                "unit": props.get("unit_of_measure"),
                "quality": props.get("approval_status"),
                "datum": None,
                "time_series_id": props.get("time_series_id"),
            })

    obs = finalize_obs(pd.DataFrame(rows))
    if derive in ("salinity", "salinity_from_conductance"):
        # PSS-78 requires the source conductance in uS/cm, so derive from the
        # unstandardized values rather than the mS/cm public product.
        obs = derive_salinity_from_usgs_conductance(obs)
    else:
        obs = standardize_usgs_units(obs)
    return obs, stations


def _normalize_unit(unit: Any) -> str:
    return (
        str(unit)
        .strip()
        .lower()
        .replace("µ", "u")
        .replace("μ", "u")
        .replace("°", "deg")
    )


def standardize_usgs_units(obs: pd.DataFrame) -> pd.DataFrame:
    """Convert common USGS continuous parameters to stable public units."""
    if obs.empty:
        return obs
    out = obs.copy()
    out["value"] = pd.to_numeric(out["value"], errors="coerce").astype(float)
    out["source_unit"] = out["unit"]
    out["unit_conversion"] = "source_unit_retained"

    for parameter, variable in USGS_STANDARD_VARIABLES.items():
        parameter_mask = out["parameter"].astype(str) == parameter
        if not parameter_mask.any():
            continue
        out.loc[parameter_mask, "variable"] = variable
        conversions = USGS_STANDARD_UNITS[parameter]
        source_units = out.loc[parameter_mask, "source_unit"].map(_normalize_unit)
        unsupported = sorted(set(source_units) - set(conversions))
        if unsupported:
            raise ValueError(
                f"Unsupported USGS unit(s) for parameter {parameter}: {unsupported}"
            )
        for source_unit, (target_unit, scale, offset) in conversions.items():
            mask = parameter_mask & (
                out["source_unit"].map(_normalize_unit) == source_unit
            )
            if not mask.any():
                continue
            out.loc[mask, "value"] = out.loc[mask, "value"] * scale + offset
            out.loc[mask, "unit"] = target_unit
            out.loc[mask, "unit_conversion"] = (
                f"value * {scale:.15g} + {offset:.15g}"
            )
    return out


def derive_salinity_from_usgs_conductance(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Convert USGS specific conductance parameter 00095 to practical salinity.

    USGS parameter 00095 is specific conductance normalized to 25 C, typically
    reported in uS/cm. The local conversion uses the PSS-78 salinometer
    polynomial with conductivity in mS/cm, temperature fixed at 25 C, and
    pressure fixed at 0 dbar. Those assumptions are recorded in provenance
    columns rather than hidden.
    """
    if obs.empty:
        return obs
    out = obs.copy()
    out["source_unit"] = out["unit"]
    unit = out["unit"].dropna().astype(str).str.lower().unique().tolist()
    if unit and not any(u in ("us/cm", "µs/cm", "umho/cm", "microsiemens per centimeter") for u in unit):
        print(f"Warning: deriving salinity from conductance with unexpected units: {unit}")

    conductance_ms_cm = pd.to_numeric(out["value"], errors="coerce") * 0.001
    out["value"] = practical_salinity_from_conductance_25c(conductance_ms_cm)
    out["variable"] = "salinity"
    out["parameter"] = "00095"
    out["unit"] = "PSU"
    out["unit_conversion"] = "PSS-78 at 25 degC and 0 dbar"
    out["quality"] = (
        out["quality"].fillna("").astype(str)
        + ";derived_from_specific_conductance_25C_p0"
    ).str.strip(";")
    out["derived_from_parameter"] = "00095"
    out["derived_assumption_temperature_C"] = 25.0
    out["derived_assumption_pressure_dbar"] = 0.0
    return out


def practical_salinity_from_conductance_25c(conductance_ms_cm: pd.Series) -> Any:
    """
    Convert conductivity in mS/cm to practical salinity at 25 C and 0 dbar.

    This is a compact PSS-78 implementation for the specific conductance use
    case: conductivity normalized to 25 C and pressure assumed 0 dbar.
    """
    values = pd.to_numeric(conductance_ms_cm, errors="coerce")

    t = 25.0
    c3515 = 42.9140
    rt35 = (
        0.6766097
        + 0.0200564 * t
        + 0.0001104259 * t**2
        - 0.00000069698 * t**3
        + 0.0000000010031 * t**4
    )
    rt = (values / c3515) / rt35
    rt = rt.clip(lower=0.0)
    rtx = rt.pow(0.5)

    a0, a1, a2, a3, a4, a5 = 0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081
    b0, b1, b2, b3, b4, b5 = 0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144
    k = 0.0162

    sp = (
        a0
        + (a1 + (a2 + (a3 + (a4 + a5 * rtx) * rtx) * rtx) * rtx) * rtx
        + ((t - 15.0) / (1.0 + k * (t - 15.0)))
        * (b0 + (b1 + (b2 + (b3 + (b4 + b5 * rtx) * rtx) * rtx) * rtx) * rtx)
    )
    return sp.where(values.notna())


def _coerce_coops_station_payload(payload: dict[str, Any]) -> dict[str, Any]:
    stations = payload.get("stations")
    if isinstance(stations, list) and stations:
        return stations[0]
    if isinstance(stations, dict):
        return stations
    return payload


def get_coops_station_info(
    station_ids: list[str],
    session: requests.Session | None = None,
) -> pd.DataFrame:
    session = session or make_session()
    rows = []
    for station_id in station_ids:
        payload = request_json(session, COOPS_STATION_URL.format(station=station_id))
        station = _coerce_coops_station_payload(payload)
        lon = station.get("lng", station.get("lon", station.get("longitude")))
        lat = station.get("lat", station.get("latitude"))
        rows.append({
            "source": "coops",
            "station_id": str(station.get("id", station_id)),
            "station_name": station.get("name"),
            "lon": pd.to_numeric(lon, errors="coerce"),
            "lat": pd.to_numeric(lat, errors="coerce"),
            "vertical_datum": None,
            "status": "ok",
        })
    return pd.DataFrame(rows)


def download_coops(
    station_ids: list[str],
    start: Any,
    end: Any,
    product: str = "water_level",
    datum: str = "NAVD",
    units: str = "metric",
    session: requests.Session | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    session = session or make_session()
    stations = get_coops_station_info(station_ids, session=session)
    station_map = stations.set_index("station_id").to_dict("index") if not stations.empty else {}
    rows = []

    for station_id in station_ids:
        params = {
            "begin_date": coops_date(start),
            "end_date": coops_date(end),
            "station": station_id,
            "product": product,
            "datum": datum,
            "time_zone": "gmt",
            "units": units,
            "application": "schism_py_pre_post",
            "format": "json",
        }
        payload = request_json(session, COOPS_DATA_URL, params=params)
        if "error" in payload:
            raise RuntimeError(f"CO-OPS {station_id}: {payload['error'].get('message')}")
        info = station_map.get(str(station_id), {})
        for item in payload.get("data", []):
            rows.append({
                "source": "coops",
                "station_id": str(station_id),
                "station_name": info.get("station_name"),
                "lon": info.get("lon"),
                "lat": info.get("lat"),
                "time": item.get("t"),
                "variable": product,
                "parameter": product,
                "value": item.get("v"),
                "unit": units,
                "quality": item.get("q"),
                "datum": datum,
                "sigma": item.get("s"),
                "flags": item.get("f"),
            })

    obs = finalize_obs(pd.DataFrame(rows))
    start_ts = parse_time(start)
    end_ts = parse_time(end)
    if start_ts is not None and end_ts is not None and not obs.empty:
        obs = obs[(obs["time"] >= start_ts) & (obs["time"] <= end_ts)]
    return obs, stations


def parse_ndbc_station_page(station_id: str, text: str) -> dict[str, Any]:
    name = None
    title = re.search(r"<title>(.*?)</title>", text, flags=re.I | re.S)
    if title:
        name = re.sub(r"\s+", " ", title.group(1)).strip()
    match = re.search(
        r"(-?\d+(?:\.\d+)?)\s+[NS]\s+(-?\d+(?:\.\d+)?)\s+[EW]",
        text,
        flags=re.I,
    )
    lat = lon = None
    if match:
        lat = float(match.group(1))
        lon = float(match.group(2))
        if re.search(rf"{match.group(1)}\s+S", text, flags=re.I):
            lat = -lat
        if re.search(rf"{match.group(2)}\s+W", text, flags=re.I):
            lon = -lon
    return {
        "source": "ndbc",
        "station_id": station_id,
        "station_name": name,
        "lon": lon,
        "lat": lat,
        "vertical_datum": None,
        "status": "ok" if lon is not None and lat is not None else "metadata_partial",
    }


def get_ndbc_station_info(
    station_ids: list[str],
    session: requests.Session | None = None,
) -> pd.DataFrame:
    session = session or make_session()
    rows = []
    for station_id in station_ids:
        response = session.get(NDBC_STATION_PAGE_URL.format(station=station_id), timeout=60)
        response.raise_for_status()
        rows.append(parse_ndbc_station_page(station_id, response.text))
    return pd.DataFrame(rows)


def read_ndbc_stdmet_text(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    if not lines:
        return pd.DataFrame()
    header = lines[0].lstrip("#").split()
    data = "\n".join(line for line in lines[2:] if line.strip() and not line.startswith("#"))
    if not data:
        return pd.DataFrame(columns=header)
    return pd.read_csv(io.StringIO(data), sep=r"\s+", names=header, na_values=["MM", "999", "9999"])


def download_ndbc(
    station_ids: list[str],
    start: Any | None = None,
    end: Any | None = None,
    variable: str = "WVHT",
    session: requests.Session | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    session = session or make_session()
    stations = get_ndbc_station_info(station_ids, session=session)
    station_map = stations.set_index("station_id").to_dict("index") if not stations.empty else {}
    start_ts = parse_time(start)
    end_ts = parse_time(end)
    variable = variable.upper()
    rows = []

    for station_id in station_ids:
        if start_ts is None or end_ts is None:
            response = session.get(NDBC_REALTIME_URL.format(station=station_id), timeout=60)
            response.raise_for_status()
            df = read_ndbc_stdmet_text(response.text)
        else:
            frames = []
            for year in range(start_ts.year, end_ts.year + 1):
                url = NDBC_HISTORICAL_URL.format(station=station_id, year=year)
                response = session.get(url, timeout=120)
                if response.status_code == 404:
                    continue
                response.raise_for_status()
                text = gzip.decompress(response.content).decode("utf-8", errors="replace")
                frames.append(read_ndbc_stdmet_text(text))
            df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

        if df.empty or variable not in df.columns:
            continue
        year_col = "#YY" if "#YY" in df.columns else "YY"
        time = pd.to_datetime(
            dict(
                year=pd.to_numeric(df[year_col], errors="coerce"),
                month=pd.to_numeric(df["MM"], errors="coerce"),
                day=pd.to_numeric(df["DD"], errors="coerce"),
                hour=pd.to_numeric(df["hh"], errors="coerce"),
                minute=pd.to_numeric(df.get("mm", 0), errors="coerce"),
            ),
            utc=True,
            errors="coerce",
        )
        info = station_map.get(station_id, {})
        for t, value in zip(time, pd.to_numeric(df[variable], errors="coerce")):
            if pd.isna(t) or pd.isna(value):
                continue
            if start_ts is not None and t < start_ts:
                continue
            if end_ts is not None and t > end_ts:
                continue
            rows.append({
                "source": "ndbc",
                "station_id": station_id,
                "station_name": info.get("station_name"),
                "lon": info.get("lon"),
                "lat": info.get("lat"),
                "time": t,
                "variable": variable,
                "parameter": variable,
                "value": value,
                "unit": None,
                "quality": None,
                "datum": None,
            })

    return finalize_obs(pd.DataFrame(rows)), stations


def finalize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "source", "station_id", "station_name", "lon", "lat", "time",
        "variable", "parameter", "value", "unit", "quality", "datum",
    ]
    for col in columns:
        if col not in obs.columns:
            obs[col] = pd.NA
    if not obs.empty:
        obs["time"] = pd.to_datetime(obs["time"], utc=True, errors="coerce")
        obs["value"] = pd.to_numeric(obs["value"], errors="coerce")
        obs["lon"] = pd.to_numeric(obs["lon"], errors="coerce")
        obs["lat"] = pd.to_numeric(obs["lat"], errors="coerce")
        obs = obs.dropna(subset=["time", "value"])
        obs = obs.sort_values(["source", "station_id", "time"]).reset_index(drop=True)
    return obs


def download_obs(
    request: ObsRequest,
    session: requests.Session | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = request.source.lower()
    if source == "coops":
        return download_coops(
            request.station_ids,
            start=request.start,
            end=request.end,
            product=request.variable or request.parameter or "water_level",
            datum=request.datum or "NAVD",
            units=request.units,
            session=session,
        )
    if source == "usgs":
        return download_usgs_continuous(
            request.station_ids,
            start=request.start,
            end=request.end,
            parameter=request.parameter or request.variable or "00060",
            derive=request.derive,
            session=session,
        )
    if source == "ndbc":
        return download_ndbc(
            request.station_ids,
            start=request.start,
            end=request.end,
            variable=request.variable or request.parameter or "WVHT",
            session=session,
        )
    raise ValueError(f"Unsupported source: {request.source}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download observations from CO-OPS, USGS, or NDBC.")
    parser.add_argument("source", choices=["coops", "usgs", "ndbc"])
    parser.add_argument("station_ids", nargs="+")
    parser.add_argument("--start")
    parser.add_argument("--end")
    parser.add_argument("--variable")
    parser.add_argument("--parameter")
    parser.add_argument("--derive", choices=["salinity", "salinity_from_conductance"])
    parser.add_argument("--datum")
    parser.add_argument("--units", default="metric")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--station-output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    obs, stations = download_obs(
        ObsRequest(
            source=args.source,
            station_ids=args.station_ids,
            start=args.start,
            end=args.end,
            variable=args.variable,
            parameter=args.parameter,
            derive=args.derive,
            datum=args.datum,
            units=args.units,
        )
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        obs.to_csv(args.output, index=False)
    else:
        print(obs.head().to_string(index=False))
        print(f"\n{len(obs)} observations")

    if args.station_output:
        args.station_output.parent.mkdir(parents=True, exist_ok=True)
        stations.to_csv(args.station_output, index=False)


if __name__ == "__main__":
    main()
