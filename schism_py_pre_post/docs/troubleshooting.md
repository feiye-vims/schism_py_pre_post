# Troubleshooting

## A USGS download appears to hang

Station discovery and API retries may be silent. HTTP 429 rate-limit responses
can trigger repeated waits. Use a USGS API key and keep it out of commands,
logs, and version control:

```bash
USGS_API_KEY="$(< /path/to/private_key_file)" python -m Download.download_usgs_polygon ...
```

If an API key is pasted into chat, a ticket, or a committed file, revoke it and
generate a replacement.

## The domain file is rejected

`download_usgs_polygon.py` requires Polygon or MultiPolygon geometry. A line
export of a grid boundary is not a polygon, even if it visually outlines the
domain. Inspect the geometry before downloading:

```bash
ogrinfo -so -al path/to/domain.gpkg
```

The reported geometry must be `Polygon` or `Multi Polygon`, and the CRS must be
defined.

## The output directory belongs to a different request

The downloader records its request in `manifest.json`. Use a new directory when
changing the polygon, date interval, parameters, derivation, or file format.
This protection prevents incompatible station parts from being combined.

## Model and BP station counts differ

The comparison script requires exactly one model column for every BP row:

```bash
awk 'NR==2 {print $1}' stations.bp
awk 'NR==1 {print NF-1; exit}' model.dat
```

Re-extract the model data with the same BP file used by the observation
workflow. Do not reorder the BP after extraction.

## MPI extraction is OOM-killed

Large SCHISM stacks can exhaust memory when too many MPI ranks read files at
the same time. Reduce the rank count and rerun:

```bash
mpiexec -n 16 python -m Extract.extract_schism_pylib_parallel ...
```

For JSON-driven plotting workflows, reduce `extract.nproc` in the run block.
Starting at half of the available cores is often more reliable than using every
core on the node.

## Matplotlib cannot write its cache

On a restricted compute node, set writable cache directories:

```bash
export MPLCONFIGDIR=/tmp/matplotlib-$USER
export XDG_CACHE_HOME=/tmp/cache-$USER
```

Then rerun the plotting command.

## A station has no observation curve

Check whether its part is empty and whether it passed the variable filter. Empty
parts can occur when metadata discovery finds a time series but the continuous
endpoint returns no values for the requested interval. Review
`selected_series.csv`, `failures.csv`, and the relevant file under
`observations/`.
