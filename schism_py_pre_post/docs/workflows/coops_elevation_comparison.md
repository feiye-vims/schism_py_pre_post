# Compare CO-OPS elevation with SCHISM runs

## Purpose

Use `Plot/plot_coops_elev.py` to compare one or more SCHISM elevation runs
against NOAA CO-OPS water-level stations. The first run in the JSON `runs`
list is the main run. It is plotted by itself as `ts_*.png`; when additional
runs are present, all runs are overlaid as `comp_*.png`.

The script also writes per-station and mean statistics. If a model time-series
file is missing, it uses `Extract/extract_schism_output.py` to extract
from raw SCHISM stacks after showing the exact MPI command and asking for
confirmation. Other plotting scripts can import the same helper.

## Prerequisites

- Python environment with the packages needed by `plot_elev.py` and
  `download_coops_elev.py`.
- Access to CO-OPS observations, or existing cached observations.
- A BP/station file whose station IDs match the CO-OPS IDs to plot.
- Existing model station time-series files, or raw SCHISM stacks under each
  run's `output_dir/outputs/`.
- For automatic extraction: `mpiexec`, `mpi4py`, SCHISM `pylib`, NumPy, and
  PyYAML. See [Extract SCHISM output at BP points with MPI](extract_schism_output_parallel.md).

## Configure a comparison

Use one JSON file per coherent comparison group. Keep shared event metadata at
the event level and run-specific paths under `runs`.

Example:

```json
{
  "All": {
    "box": {"W": -100, "E": -60, "S": 8, "N": 48}
  },
  "events": {
    "2018_hindcast": {
      "plot_start_day_str": "2017-12-01 00:00:00",
      "plot_end_day_str": "2018-04-01 00:00:00",
      "station_bp_file": "/sciclone/schism10/feiye/STOFS3D-v7.3/R23a/station.in",
      "datum_shift_file": "/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/navd2xgeoid_shift.txt",
      "default_datum": "NAVD",
      "runs": [
        {
          "name": "2018_v7p3_2D",
          "model_start_day_str": "2017-12-01 00:00:00",
          "elev_out_file": "/sciclone/schism10/feiye/STOFS3D-v7.3/O23e/staout_1",
          "output_dir": "/sciclone/schism10/feiye/STOFS3D-v7.3/O23e/",
          "shift": 0.0,
          "line_style": "k"
        },
        {
          "name": "R15b4_v7",
          "model_start_day_str": "2017-12-01 00:00:00",
          "station_bp_file": "/sciclone/schism10/feiye/STOFS3D-v7.3/BPfiles/stofs3d_atl_202503.bp",
          "elev_out_file": "/sciclone/schism10/feiye/STOFS3D-v8/O15b4_v7/elevation.stofs3d_atl_202503.dat",
          "output_dir": "/sciclone/schism10/feiye/STOFS3D-v8/O15b4_v7/",
          "shift": 0.0,
          "line_style": "--g",
          "extract": {
            "var_name": "elevation",
            "nproc": 4
          }
        }
      ]
    }
  }
}
```

Important fields:

- `station_bp_file`: event-level station file used for station grouping and
  observations. A run-level `station_bp_file` overrides it only for that run's
  model file.
- `elev_out_file`: model time-series file consumed by `plot_elev.py`.
- `output_dir`: SCHISM run/output directory. If `elev_out_file` is missing,
  the script expects raw stacks under `output_dir/outputs/`.
- `extract`: optional extraction overrides. Do not set `run_dir`,
  `output_file`, `start_stack`, or `end_stack`; these are inferred.

## Run the comparison

From the repository root:

```bash
MPLCONFIGDIR=/tmp python -m Plot.plot_coops_elev \
  Plot/coops_elev_2018_v7_runs.json \
  2018_hindcast \
  --region Full_domain \
  --output-dir .
```

`MPLCONFIGDIR=/tmp` avoids Matplotlib cache-permission warnings on systems
where the default home config directory is not writable.

Useful options:

- `--no-scatter`: skip map/scatter plots, useful when Basemap is unavailable.
- `--no-extract-missing`: fail instead of extracting missing model files.
- `--extract-mpi-np N`: default MPI rank count for extraction.
- `--extract-script PATH`: use a non-default extraction script path.
- `--extract-config-dir PATH`: directory for generated extraction configs.

Choose `--extract-mpi-np` based on node memory as well as CPU count. A practical
starting point on memory-constrained nodes is half of `nproc`; if MPI reports an
OOM kill, reduce the rank count and rerun. Successful extraction publishes the
destination file only after the write completes, so a killed run should not
leave a partial target file.

## Missing model extraction

When `elev_out_file` is missing or empty, the script validates
`output_dir/outputs/` through `Extract/extract_schism_output.py` before
proposing extraction. Exactly one of these layouts must be present:

```text
schout_1.nc, schout_2.nc, ...
out2d_1.nc, out2d_2.nc, ...
```

The stack numbers must be consecutive. File sizes are checked before
extraction. Minor size differences are expected, so the helper compares the
final stack against the mean and standard deviation of the previous stacks. If
the final stack is a strong low-size outlier, only that trailing stack is
excluded from extraction. When the previous stacks are effectively identical in
size, the helper falls back to a 1% low-size tolerance.

Before running MPI, the script prints the raw-output directory, stack type,
included stacks, excluded stacks, generated config path, and exact command:

```text
Stacks to combine: 1-35
Stacks excluded: 36
Extraction command:
mpiexec -n 4 python -m Extract.extract_schism_pylib_parallel --config ... --run R15b4_v7
Proceed with extraction? [y/N]:
```

After confirmation, the script checks `elev_out_file` again. If the file
appeared while the prompt was waiting, extraction is skipped and the existing
file is used.

To reuse this behavior in another plotting script, import `ensure_model_file`:

```python
from Extract.extract_schism_output import ensure_model_file

ensure_model_file(
    run_spec,
    extract_missing=True,
    extract_mpi_np=4,
    extract_config_dir="extract_configs",
)
```

`run_spec` should contain the same run-level fields used in the comparison
JSON: `name`, `elev_out_file`, `output_dir`, `station_bp_file`, and optionally
`extract`.

## Outputs

Typical outputs include:

- `ts_<event>_<group>_<datum>.png`: observations and the main run.
- `comp_<event>_<group>_<datum>.png`: observations and all configured runs.
- `stats_<run>_<group>.txt` and `.csv`: per-station stats for a station group.
- `stats_<run>_<event>_<region>.txt` and `.csv`: overall stats.
- `mean_stats_<run>_<event>_<region>.txt`: group means and overall mean.
- `datum_info_<event>_<region>.txt`: count of stations using NAVD, MSL, or no data.
- `extract_configs/extract_<run>.yaml`: generated extraction config when
  missing-model extraction is requested.

## Validation

Validate the JSON before running:

```bash
python -m json.tool Plot/coops_elev_2018_v7_runs.json >/tmp/coops_elev.json
```

Check that run labels and line styles are read as expected:

```bash
MPLCONFIGDIR=/tmp python - <<'PY'
from Plot.plot_coops_elev import _expand_run_specs
_, runs = _expand_run_specs("Plot/coops_elev_2018_v7_runs.json", "2018_hindcast")
print([run["name"] for run in runs])
print([run["line_style"] for run in runs])
PY
```

After plotting, confirm that expected stats and plots exist:

```bash
ls ts_2018_hindcast_*png comp_2018_hindcast_*png stats_*2018_hindcast*txt
```

## Assumptions and limitations

- CO-OPS observation times are handled as UTC.
- `model_start_day_str` must match the model time origin in `elev_out_file`.
- The station IDs in each model file must include the stations being plotted.
- Run-level station files are allowed, but their station IDs must still match
  the comparison station IDs.
- Missing-model extraction supports scalar variables accepted by
  `pylib.read_schism_output`.
- Raw-output auto-detection accepts only `schout_*.nc` or `out2d_*.nc`, not a
  mixture of both.
- Only a smaller trailing stack that is a strong file-size outlier is
  auto-excluded. Any other missing or non-consecutive stack condition is an
  error.

See [general troubleshooting](../troubleshooting.md) for common failures.
