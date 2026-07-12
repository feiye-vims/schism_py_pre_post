# Extract SCHISM output at BP points with MPI

## Purpose

Use `extract_schism_pylib_parallel.py` to extract a scalar SCHISM variable at
the points in a BP file. The script divides an inclusive range of SCHISM output
stacks among MPI ranks, gathers the results, and writes one whitespace-delimited
time-series file.

This workflow is useful for producing station files consumed by plotting and
model/observation comparison scripts. Plotting workflows can also call
`Extract/extract_schism_output.py`, which validates raw stack layout,
generates the YAML configuration, prints the `mpiexec` command, and asks for
confirmation before running this MPI extractor.

## Prerequisites

- A Python environment containing `mpi4py`, NumPy, PyYAML, and SCHISM `pylib`.
- An MPI implementation compatible with the installed `mpi4py` package.
- A completed SCHISM run with readable output stacks.
- A BP file containing the extraction points in the required output-column
  order.
- A YAML configuration describing one or more extraction cases.

The extraction module used in this example is:

```text
Extract.extract_schism_pylib_parallel
```

## Configure extraction cases

By default, the script reads:

```text
/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/extract.yaml
```

The configuration may be a list of cases or a mapping with a `runs` list. Each
case needs a unique `name` when it will be selected with `--run`:

```yaml
runs:
  - name: R19i1
    run_dir: /sciclone/schism10/feiye/STOFS3D-v8/R19i1/
    var_name: salinity
    output_file: /sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.usgs_sal_la.dat
    bpfile: /sciclone/schism10/feiye/STOFS3D-v8/BPfiles/usgs_sal_la.bp
    start_stack: 1
    end_stack: 35
```

The stack interval is inclusive. `var_name` must be a scalar variable supported
by `pylib.read_schism_output`. The output parent directory is created if it does
not exist.

## Run one configured case

Run the `R19i1` case with four MPI ranks:

```bash
mpiexec -n 4 python -m Extract.extract_schism_pylib_parallel \
  --run R19i1
```

The script splits stacks `1` through `35` into four contiguous rank-local
ranges. Rank 0 gathers the extracted arrays and writes the configured output
file after all ranks finish successfully.

To use a configuration other than the script default, pass it explicitly:

```bash
mpiexec -n 4 python -m Extract.extract_schism_pylib_parallel \
  --config /path/to/extract.yaml \
  --run R19i1
```

Specify `--run` more than once to process selected cases in the requested
order. Omit `--run` to process every case in the YAML file:

```bash
mpiexec -n 4 python -m Extract.extract_schism_pylib_parallel \
  --config /path/to/extract.yaml \
  --run R19i1 \
  --run R29i1
```

Useful safeguards and formatting options are:

- `--no-overwrite`: fail if the destination file already exists. Without this
  option, a successful extraction replaces the destination atomically.
- `--fmt FORMAT`: override the NumPy `savetxt` format; the default is `%.4f`.
- `--single-case`: ignore YAML and use the constants defined near the top of
  the script. This mode cannot be combined with `--run`.

Choose the MPI rank count based on memory as well as CPU count. More ranks read
more stack files concurrently, which can increase memory pressure. A practical
first retry after an OOM kill is half of the available cores, then halve again
if needed.

## Output format

The output is a numeric text table without a header:

```text
elapsed_time point_1 point_2 ... point_N
```

The first column is the time returned by `pylib`; each remaining column is the
scalar value at one BP point, in BP-file order. The script writes to a temporary
file in the destination directory and publishes it only after the write
completes.

For the `R19i1` example, the output is:

```text
/sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.usgs_sal_la.dat
```

## Validate the extraction

Confirm that the file exists and is nonempty:

```bash
test -s /sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.usgs_sal_la.dat
```

Confirm that its value-column count equals the BP station count:

```bash
awk 'NR==2 {print "BP stations:", $1}' \
  /sciclone/schism10/feiye/STOFS3D-v8/BPfiles/usgs_sal_la.bp

awk 'NR==1 {print "Model stations:", NF-1; exit}' \
  /sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.usgs_sal_la.dat
```

Inspect the first and last output records to check time coverage:

```bash
head -n 1 /sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.usgs_sal_la.dat
tail -n 1 /sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.usgs_sal_la.dat
```

## Next step

To compare extracted salinity with observations, follow
[Compare USGS station variables with SCHISM](usgs_station_variable_comparison.md).
That workflow owns the observation download and comparison commands for gage
height, salinity, temperature, and similar station variables.

## Assumptions and limitations

- Only scalar extracted data are accepted. Vector or multi-level variables do
  not match the expected two-dimensional point-by-time layout.
- All MPI ranks must use the same environment and be able to read the run and
  BP paths and write the output directory.
- More MPI ranks than stacks is allowed, but excess ranks do no extraction and
  provide no speed benefit.
- The output contains no header or BP metadata. Preserve the YAML configuration
  and BP file with the result.
- Correct model start time and time units must be supplied to downstream tools;
  they are not recorded in the extracted file.

## Troubleshooting

- If `--run` reports that a name is missing, check the spelling and ensure the
  YAML case has a unique `name`.
- If validation fails before extraction, verify `run_dir`, `bpfile`, and the
  inclusive stack range.
- If one rank fails, the script reports errors collected from all ranks and
  does not publish a partial output file.
- If MPI reports an OOM kill, reduce `mpiexec -n` and rerun. The destination
  file is written through a temporary file and should appear only after a
  successful run.
- If MPI fails during import or startup, ensure `mpiexec`, Python, and `mpi4py`
  come from compatible installations.

See [general troubleshooting](../troubleshooting.md) for additional failures.
