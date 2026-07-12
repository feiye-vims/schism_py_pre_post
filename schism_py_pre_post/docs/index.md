# Documentation

This directory contains practical, reproducible workflows for this repository.
The emphasis is on commands that can be rerun, expected file formats, validation
checks, and common failure modes.

## Workflows

- [Compare CO-OPS elevation with SCHISM runs](workflows/coops_elevation_comparison.md)
- [Extract SCHISM output at BP points with MPI](workflows/extract_schism_output_parallel.md)
- [Compare USGS station variables with SCHISM](workflows/usgs_station_variable_comparison.md)
- [Template for a new workflow](workflows/template.md)

## Reference

- [Troubleshooting](troubleshooting.md)

## Documentation conventions

Use one Markdown file per workflow. A workflow should include:

1. Purpose and scope.
2. Required software and input files.
3. Commands with all important options shown.
4. Output files and their formats.
5. Quick validation checks.
6. Known assumptions and limitations.
7. Links to relevant scripts.

Use paths relative to the repository when referring to source code. Absolute
paths are acceptable in site-specific command examples, but label them as
examples so they are not mistaken for portable defaults.
