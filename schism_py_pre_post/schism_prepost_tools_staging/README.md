# SCHISM Pre/Post Tools Staging

This folder is a symlink-based staging area for a curated, shareable subset of
the main repository. It is not yet a standalone package.

## Included Scripts

- `Download/download_obs.py`
- `Download/download_usgs_polygon.py`
- `Extract/extract_schism_output.py`
- `Extract/extract_schism_pylib_parallel.py`
- `Plot/compare_usgs_schism.py`
- `Plot/plot_coops_elev.py`
- `Stats/model_obs_io.py`
- `Stats/model_obs_statistics.py`

## Notes

- The files here are symlinks to the main repository sources, so edits should
  be made in the canonical source locations.
- The folder names intentionally mirror the main repository to avoid changing
  imports during staging.
- Run entry points as modules from this folder or the main repository root, for
  example `python -m Plot.compare_usgs_schism` or
  `python -m Extract.extract_schism_output`.
- Before publishing or making a pip package, export this folder with symlinks
  dereferenced, for example:

```bash
rsync -aL schism_prepost_tools_staging/ /tmp/schism_prepost_tools_release/
```
