#!/usr/bin/env python

import os
os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg", force=True)

import numpy as np
import xarray as xr
from mpi4py import MPI

from pylib import read_schism_output, read

# --------------- inputs ----------------
RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R39b/'
output_files = [
    '/sciclone/schism10/feiye/STOFS3D-v8/O39b/TS.transect.mississippi2',
]
bpfile_list = [
    '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/mississippi2.bp'
]

start_stack = 1
end_stack = 80

EXTRACT_VARS = ['salinity', 'zCoordinates']
# EXTRACT_VARS = ['viscosity', 'diffusivity', 'verticalVelocity', 'horizontalVelX', 'horizontalVelY']
# ---------------------------------------


def extract_transect(stack_range=None, bpfile=None):
    """
    Read SCHISM transect for the given stack_range.
    """
    data = read_schism_output(
        run=RUN_DIR,
        varname=EXTRACT_VARS,
        xyz=bpfile,
        stacks=stack_range,
        fmt=1,  # transect xy and all levels
    )
    return data


def result_to_piece_dict(result, stack_start, stack_end):
    """
    Convert pylib result object into a plain dict that MPI can gather reliably.
    """
    # INFO is used only to infer variable names in the original script;
    # we keep it to stay close to your original logic.
    info = list(getattr(result, "INFO", []))

    piece = {
        "stack_start": int(stack_start),
        "stack_end": int(stack_end),
        "INFO": info,
        "time": np.asarray(result.time),
    }

    # Known vars for transect mode in this script
    for v in EXTRACT_VARS:
        if hasattr(result, v):
            piece[v] = np.asarray(getattr(result, v))
        else:
            raise AttributeError(f"Result missing attribute '{v}'. Available: {dir(result)}")

    return piece


def combine_piece_dicts(pieces):
    """
    Combine gathered pieces (dicts) into a single xarray Dataset.
    This mirrors your original combine_results() but avoids pickling pylib objects.
    """
    # Ensure stack order is preserved
    pieces = sorted(pieces, key=lambda d: d["stack_start"])

    # Variable name discovery (fallback to known list)
    if pieces and pieces[0].get("INFO"):
        variable_names = [s.split(":")[0].strip() for s in pieces[0]["INFO"]]
    else:
        variable_names = ["time"] + EXTRACT_VARS

    ds = xr.Dataset()

    for var_name in variable_names:
        if var_name == "time":
            t = np.concatenate([p["time"] for p in pieces], axis=0)
            ds["time"] = (("time",), t)
        elif var_name in ["elevation", "elev"]:
            var_combined = np.concatenate([p[var_name] for p in pieces], axis=1)
            data = np.transpose(var_combined)  # -> (nt, npoints)
            ds[var_name] = (("time", "npoints"), data)
        elif var_name in ["salinity", "temperature", "zCoordinates", "viscosity", "diffusivity", "verticalVelocity", "horizontalVelX", "horizontalVelY"]:
            # In pylib, these are typically shaped like (npoints, nt, nz) per chunk,
            var_combined = np.concatenate([p[var_name] for p in pieces], axis=1)  # concat along nt-dim
            data = np.transpose(var_combined, (1, 0, 2))  # -> (nt, npoints, nz)
            ds[var_name] = (("time", "npoints", "nz"), data)
        else:
            raise ValueError(f"Unexpected variable '{var_name}'")

    return ds


def mpi_extract():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Build stack ranges and split across ranks (keeps original "array_split" behavior)
    all_stacks = np.arange(start_stack, end_stack + 1)
    stack_ranges = np.array_split(all_stacks, size)
    my_range = stack_ranges[rank]

    # Each rank reads its chunk (or returns None if empty)
    my_piece = None
    if my_range.size > 0:
        if rank == 0:
            print(f"[rank {rank}/{size}] total stacks={len(all_stacks)}; splitting into {size} ranks")

        print(f"[rank {rank}] Processing stacks {int(my_range[0])} to {int(my_range[-1])} (n={len(my_range)})")

        result = extract_transect(stack_range=my_range, bpfile=bpfile_list[0])
        my_piece = result_to_piece_dict(result, stack_start=my_range[0], stack_end=my_range[-1])
    else:
        print(f"[rank {rank}] No stacks assigned.")

    # Gather pieces to root
    pieces = comm.gather(my_piece, root=0)

    # Root combines and writes
    if rank == 0:
        pieces = [p for p in pieces if p is not None]

        if not pieces:
            raise RuntimeError("No results gathered from any rank.")

        results_ds = combine_piece_dicts(pieces)

        # attach lon/lat from bpfile
        bp = read(bpfile_list[0])
        results_ds["lon"] = (("npoints",), bp.x)
        results_ds["lat"] = (("npoints",), bp.y)

        out_nc = output_files[0] + ".nc"
        results_ds.to_netcdf(out_nc)
        print(f"[rank 0] Wrote: {out_nc}")


if __name__ == "__main__":
    mpi_extract()
    # Make sure everyone finishes cleanly before exit
    MPI.COMM_WORLD.Barrier()
    if MPI.COMM_WORLD.Get_rank() == 0:
        print("done!")
