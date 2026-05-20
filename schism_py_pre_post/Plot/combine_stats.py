import os
import re
from glob import glob
from typing import List, Dict, Optional

import pandas as pd

# ---------------- user inputs ----------------
input_glob = "/sciclone/schism10/feiye/STOFS3D-v7.3/Total/stats_Total_*.txt"   # <-- change me (or "*.txt" in current dir)
out_path = "/sciclone/schism10/feiye/STOFS3D-v7.3/Total/combined_means_stats"  # <-- change me
sort_regions = True                  # optional
neglected_stations = [
    "8632837",  # Rappahannock Light    
]
# --------------------------------------------

COLUMNS = [
    "station_name",
    "station_id",
    "station_lon",
    "station_lat",
    "RMSE",
    "MAE",
    "Bias",
    "CC",
    "ubRMSE",
    "Max_Obs",
    "Max_Mod",
]
METRICS = ["RMSE", "MAE", "Bias", "CC", "ubRMSE", "Max_Obs", "Max_Mod"]
META_KEEP = ["station_lon", "station_lat"]

# First token that is >= 4 digits is station_id; everything before is station_name
ID_TOKEN_RE = re.compile(r"\s+(\d{4,})\s+")


def infer_region_name(path: str) -> str:
    """Use filename stem as region label."""
    return os.path.splitext(os.path.basename(path))[0]


def to_float(tok: str) -> float:
    t = tok.strip()
    if t.lower() == "nan":
        return float("nan")
    return float(t)


def parse_data_line(line: str) -> Optional[Dict[str, object]]:
    """
    Parse one line using the rule:
      - station_name: everything before the first >=4-digit token
      - station_id: that token
      - remaining columns: whitespace-separated
    Returns None for header/blank/non-data lines.
    """
    line = line.rstrip("\n")
    if not line.strip():
        return None

    # Skip header line(s)
    if line.lstrip().startswith("station_name"):
        return None

    m = ID_TOKEN_RE.search(line)
    if not m:
        return None

    station_id = m.group(1)
    station_name = line[: m.start()].strip()  # removes leading spaces
    rest = line[m.end():].strip()

    rest_tokens = rest.split()
    expected = len(COLUMNS) - 2  # excluding station_name and station_id

    if len(rest_tokens) != expected:
        raise ValueError(
            f"Bad field count while parsing.\n"
            f"File line: {line}\n"
            f"Parsed station_name={station_name!r}, station_id={station_id!r}\n"
            f"Expected {expected} trailing fields, got {len(rest_tokens)}: {rest_tokens}"
        )

    vals = [to_float(t) for t in rest_tokens]
    row: Dict[str, object] = {"station_name": station_name, "station_id": station_id}
    for k, v in zip(COLUMNS[2:], vals):
        row[k] = v
    return row


def read_station_table(path: str) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parsed = parse_data_line(line)
            if parsed is not None:
                rows.append(parsed)

    if not rows:
        raise ValueError(f"{path}: no data rows parsed (check format).")

    df = pd.DataFrame(rows, columns=COLUMNS)

    # Drop per-file Mean row(s)
    df = df[df["station_name"] != "Mean"].copy()

    # Numeric safety (already floats, but keep robust)
    for c in METRICS + META_KEEP:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Keep station_id as string (safe)
    df["station_id"] = df["station_id"].astype(str).str.strip()

    return df


def mean_row(df: pd.DataFrame, region_label: str) -> Dict[str, object]:
    out: Dict[str, object] = {
        "region": region_label,
        "nstation": int(len(df)),
    }
    out.update(df[METRICS].mean(skipna=True).to_dict())
    for c in META_KEEP:
        out[c] = df[c].mean(skipna=True)
    return out


# ---------------- main workflow ----------------
files = sorted(glob(input_glob))
if not files:
    raise FileNotFoundError(f"No files matched: {input_glob}")

per_region_rows: List[Dict[str, object]] = []
all_station_dfs: List[pd.DataFrame] = []

for path in files:
    df = read_station_table(path)
    # Skip neglected stations
    df = df[~df["station_id"].isin(neglected_stations)].copy()

    all_station_dfs.append(df)

    region = infer_region_name(path)
    per_region_rows.append(mean_row(df, region))

# pooled ALL stations (NOT mean of regional means)
df_all = pd.concat(all_station_dfs, ignore_index=True)
all_row = mean_row(df_all, "ALL")

out_df = pd.DataFrame(per_region_rows)
if sort_regions:
    out_df = out_df.sort_values("region", kind="stable")

out_df = pd.concat([out_df, pd.DataFrame([all_row])], ignore_index=True)

# Output columns: only regional means + pooled mean
cols = ["region", "nstation"] + METRICS
out_df = out_df[cols]

# Write aligned text
os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
with open(out_path, "w", encoding="utf-8") as f:
    f.write(
        out_df.to_string(
            index=False,
            justify="right",
            float_format=lambda x: f"{x:10.4f}",
        )
    )
    f.write("\n")

print(f"Wrote: {out_path}")
print(f"Regions: {len(files)} | Total pooled stations: {len(df_all)}")
