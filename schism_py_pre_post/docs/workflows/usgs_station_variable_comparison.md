# Compare USGS station variables with SCHISM

## Purpose

This workflow downloads USGS observations inside a polygon, extracts matching
SCHISM station time series, and plots model/observation comparisons. It works
for gage height, salinity, temperature, and other scalar USGS parameters
supported by `Download/download_usgs_polygon.py`.

For salinity, `salinity_best_available` applies this selection rule:

1. Use observed USGS salinity, parameter `00480`, when it exists.
2. Otherwise use specific conductance, parameter `00095`, and derive practical
   salinity.
3. If both parameters exist, retain only observed salinity.

The downloader records the choice in `observation_type` and `derivation`.

## Prerequisites

- Python with pandas, GeoPandas, PyArrow, Matplotlib, and Requests.
- A USGS API key stored in a file readable only by the user.
- A Polygon or MultiPolygon domain in a GeoPackage, Shapefile, or GeoJSON.
- SCHISM raw stacks under each run's `output_dir/outputs/`, or existing
  extracted model station files.

The commands below are a template. Replace the shell variables with the
polygon, date range, run directory, extracted station-file path, and output
locations for the case being compared.

## 1. Download Observations

The example key file should have mode `600`:

```bash
chmod 600 /sciclone/home/feiye/usgs_api_key
```

Set the common inputs:

```bash
POLYGON=/path/to/domain_polygon.gpkg
START_DATE=YYYY-MM-DD
END_DATE=YYYY-MM-DD
MODEL_START="YYYY-MM-DD 00:00:00"
OBS_ROOT=/path/to/observations
RUN_DIR=/path/to/schism_run
POST_DIR=/path/to/postprocessed_outputs
COMPARISON_ROOT=/path/to/comparison_outputs
CONFIG_JSON=/path/to/usgs_variable_comparison.json
```

Download best-available salinity:

```bash
USGS_API_KEY="$(< /sciclone/home/feiye/usgs_api_key)" \
python -m Download.download_usgs_polygon \
  "$POLYGON" \
  --start "$START_DATE" \
  --end "$END_DATE" \
  --parameter salinity_best_available \
  --output-dir "$OBS_ROOT/usgs_salinity_best_${START_DATE//-/}_${END_DATE//-/}"
```

Download water temperature:

```bash
USGS_API_KEY="$(< /sciclone/home/feiye/usgs_api_key)" \
python -m Download.download_usgs_polygon \
  "$POLYGON" \
  --start "$START_DATE" \
  --end "$END_DATE" \
  --parameter temperature \
  --output-dir "$OBS_ROOT/usgs_temperature_${START_DATE//-/}_${END_DATE//-/}"
```

The output is resumable. Rerunning the identical request skips completed
station/parameter parts. Do not reuse an output directory for a different
polygon, date range, parameter, derivation method, or output format.

Important outputs:

| Path | Description |
|---|---|
| `observations/*.parquet` | Full observation time series, one station/parameter part per file. |
| `discovered_series.csv` | All series returned during station discovery. |
| `selected_series.csv` | Series retained for download. |
| `stations.csv` | Station metadata and coordinates. |
| `mean_xyz.csv` | Station means with variable, unit, and derivation metadata. |
| `mean_xyz.bp` | SCHISM BP points in selected-series order. |
| `manifest.json` | Request definition and completion counts. |
| `failures.csv` | Station/parameter requests that failed. |

## 2. Configure Model Comparisons

Put plot and extraction settings in a JSON file. Example:

```json
{
  "comparisons": [
    {
      "name": "salinity_run_name_USGS",
      "obs_dir": "/path/to/observations/usgs_salinity_best_YYYYMMDD_YYYYMMDD",
      "bp_file": "/path/to/observations/usgs_salinity_best_YYYYMMDD_YYYYMMDD/mean_xyz.bp",
      "model_start": "YYYY-MM-DD 00:00:00",
      "model_time_unit": "days",
      "variable": "salinity",
      "ylabel": "Salinity (PSU)",
      "start": "YYYY-MM-DD",
      "end": "YYYY-MM-DD",
      "overview_layout": [6, 3],
      "overview_ylim": [0, 35],
      "output_dir": "/path/to/comparison_outputs/salinity_run_name_USGS",
      "runs": [
        {
          "name": "run_name",
          "output_dir": "/path/to/schism_run",
          "elev_out_file": "/path/to/postprocessed_outputs/salinity.usgs.dat",
          "extract": {
            "var_name": "salinity",
            "nproc": 4
          }
        }
      ]
    }
  ]
}
```

`output_dir` in a run is the SCHISM run directory containing `outputs/` and
mesh metadata such as `vgrid.in`. `elev_out_file` is the extracted station file
to plot; the historical name is kept for compatibility even when the variable
is salinity or temperature.

## 3. Extract And Plot

Run all comparisons from the JSON:

```bash
MPLCONFIGDIR=/tmp python -m Plot.compare_usgs_schism \
  --config "$CONFIG_JSON" \
  --extract-config-dir /tmp/usgs_extract_configs
```

Run only one comparison:

```bash
MPLCONFIGDIR=/tmp python -m Plot.compare_usgs_schism \
  --config "$CONFIG_JSON" \
  --comparison salinity_run_name_USGS \
  --extract-config-dir /tmp/usgs_extract_configs
```

If a model station file is missing, `compare_usgs_schism.py` calls
`Extract/extract_schism_output.py`. The helper checks the variable's
own raw stack files:

- salinity uses `salinity_*.nc`
- temperature uses `temperature_*.nc`
- elevation/gage height uses `out2d_*.nc`, with `schout_*.nc` as fallback

It validates consecutive stacks, checks the trailing stack size against the
mean and standard deviation of previous stacks, prints the exact `mpiexec`
command, and asks before running extraction.

Set each run's `extract.nproc` according to available memory as well as CPU
count. If a full-core extraction is OOM-killed, halve `nproc` and rerun. A
successful extraction writes the target station file only after the MPI run
finishes, so a failed attempt should not leave a partial target file.

## Outputs

The comparison script writes:

- One auto-scaled PNG for each station.
- Paginated overview PNGs using the requested subplot layout.
- Overview images at 3840 x 2160 pixels with a 16:9 aspect ratio.

Outputs go to each comparison's `output_dir`.

## Validation

Validate the JSON:

```bash
python -m json.tool "$CONFIG_JSON" >/tmp/usgs_variable_comparison.json
```

Confirm observation coverage:

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd

for obs_dir in [
    "/path/to/observations/usgs_salinity_best_YYYYMMDD_YYYYMMDD",
    "/path/to/observations/usgs_temperature_YYYYMMDD_YYYYMMDD",
]:
    parts = Path(obs_dir, "observations").glob("*.parquet")
    frames = [pd.read_parquet(path) for path in parts]
    obs = pd.concat([frame for frame in frames if not frame.empty], ignore_index=True)
    obs["time"] = pd.to_datetime(obs["time"], utc=True)
    print(obs_dir)
    print("rows:", len(obs))
    print("stations:", obs["station_id"].nunique())
    print("time range:", obs["time"].min(), obs["time"].max())
PY
```

## Assumptions And Limitations

- Observation timestamps are normalized to UTC.
- `model_start` must match the model time origin in the extracted model file.
- Station identity and model-column order come from `bp_file`.
- Run-level raw stack detection is variable-specific; salinity and temperature
  do not use `out2d_*.nc` for file-size checks.
- These observations do not contain sensor-depth information. The comparison
  is therefore not depth-matched in stratified water.
- Conductance-derived salinity uses the derivation recorded by the downloader.
- The workflow plots time series but does not compute skill metrics; use
  `python -m Stats.model_obs_statistics` for statistics and diagrams.

See [general troubleshooting](../troubleshooting.md) for common failures.
