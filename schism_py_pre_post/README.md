# SCHISM Python Pre/Post-Processing

Utilities and workflows for preparing SCHISM inputs, downloading observation
data, extracting model results, and comparing model output with observations.

## Main Workflows

- Compare SCHISM elevation with NOAA CO-OPS stations:
  [docs/workflows/coops_elevation_comparison.md](docs/workflows/coops_elevation_comparison.md)
- Download USGS station observations and compare gage height, salinity,
  temperature, or other scalar variables:
  [docs/workflows/usgs_station_variable_comparison.md](docs/workflows/usgs_station_variable_comparison.md)
- Extract SCHISM output at BP points with MPI:
  [docs/workflows/extract_schism_output_parallel.md](docs/workflows/extract_schism_output_parallel.md)

## Common Entry Points

- `python -m Plot.plot_coops_elev`: CO-OPS elevation comparison from a JSON
  file.
- `python -m Plot.compare_usgs_schism`: USGS station-variable comparison from a
  JSON file or legacy command-line arguments.
- `python -m Download.download_usgs_polygon`: USGS observation download inside
  a Polygon or MultiPolygon domain.
- `python -m Extract.extract_schism_output`: reusable helper that validates raw
  SCHISM stacks, generates extraction configs, prints the MPI command, and asks
  before extracting missing model station files.
- `python -m Extract.extract_schism_pylib_parallel`: MPI extractor for scalar
  SCHISM variables at BP points.
- `python -m Stats.model_obs_statistics`: model/observation metrics plus
  target and Taylor diagrams.
- `Stats/model_obs_io.py`: shared model/observation readers used by plotting
  and statistics modules.

## Documentation

- [Documentation index](docs/index.md)
- [CO-OPS elevation comparison](docs/workflows/coops_elevation_comparison.md)
- [USGS station-variable comparison](docs/workflows/usgs_station_variable_comparison.md)
- [Parallel SCHISM extraction](docs/workflows/extract_schism_output_parallel.md)
- [Workflow documentation template](docs/workflows/template.md)
- [Troubleshooting](docs/troubleshooting.md)

The files under `docs/workflows/` are the canonical instructions for regularly
used procedures. Keep machine-specific paths in clearly marked examples and
explain the expected inputs and outputs for every command.

## Notes

- Set `MPLCONFIGDIR=/tmp` or another writable directory when plotting on
  restricted compute nodes.
- JSON-driven plotting workflows can extract missing model station files from
  `$output_dir/outputs/` when raw SCHISM stacks are present.
- Choose MPI rank counts based on available memory as well as CPU count; if an
  extraction is OOM-killed, reduce `nproc` and rerun.
