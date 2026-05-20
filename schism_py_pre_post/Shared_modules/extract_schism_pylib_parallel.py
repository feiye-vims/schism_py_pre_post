#!/usr/bin/env python3
from pathlib import Path

import numpy as np
from mpi4py import MPI
import yaml

from pylib import read_schism_output


# --------------------------- top-level inputs -----------------------
# Set to None for original single-case mode
# CONFIG_FILE = None
CONFIG_FILE = None  # '/sciclone/home/feiye/spp/Shared_modules/extract.yaml'

if CONFIG_FILE is None:
    # --------------- single-case inputs ----------------
    RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v7.3/F2018/'
    VAR_NAME = 'elev'
    OUTPUT_FILE = '/sciclone/schism10/feiye/STOFS3D-v7.3/F2018/elevation.stofs3d_atl_202501.dat'
    BPFILE = '/sciclone/schism10/feiye/STOFS3D-v7.3/BPfiles/stofs3d_atl_202501.bp'

    start_stack = 2
    end_stack = 366

np_savetxt_args = {'fmt': '%.4f', 'delimiter': ' ', 'newline': '\n'}
# ---------------------------------------------------


def process_stack_range(run_dir: str, var_name: str, stack_range: np.ndarray, bpfile: str) -> np.ndarray:
    """
    Process a range of stacks (inclusive indices provided by stack_range array).
    Returns a 2D array: [time_col, data_cols...]
    """
    if stack_range.size == 0:
        return np.empty((0, 0), dtype=float)

    rank = MPI.COMM_WORLD.Get_rank()
    print(
        f"[rank {rank}] processing stacks {int(stack_range[0])}..{int(stack_range[-1])} "
        f"for run={run_dir}",
        flush=True
    )

    data = read_schism_output(
        run=run_dir,
        varname=[var_name],
        xyz=bpfile,
        stacks=stack_range
    )

    out = np.c_[data.time, getattr(data, var_name).T]
    return out


def split_stacks(start_s: int, end_s: int, size: int, rank: int) -> np.ndarray:
    """
    Split [start_s, end_s] into 'size' contiguous chunks and return this rank's chunk.
    """
    stacks = np.arange(start_s, end_s + 1, dtype=int)
    chunks = np.array_split(stacks, size)
    return chunks[rank]


def process_one_case(case: dict):
    """
    Required keys:
      - run_dir
      - var_name
      - output_file
      - bpfile
      - start_stack
      - end_stack
    Optional keys:
      - name
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    run_dir = case['run_dir']
    var_name = case['var_name']
    output_file = case['output_file']
    bpfile = case['bpfile']
    start_s = int(case['start_stack'])
    end_s = int(case['end_stack'])
    case_name = case.get('name', 'unnamed_case')

    if rank == 0:
        print(f"\n=== Processing case: {case_name} ===", flush=True)
        print(f"RUN_DIR    : {run_dir}", flush=True)
        print(f"VAR_NAME   : {var_name}", flush=True)
        print(f"BPFILE     : {bpfile}", flush=True)
        print(f"STACKS     : {start_s}..{end_s}", flush=True)
        print(f"OUTPUT     : {output_file}", flush=True)

    comm.Barrier()

    local_stacks = split_stacks(start_s, end_s, size, rank)
    local_result = process_stack_range(run_dir, var_name, local_stacks, bpfile)

    gathered = comm.gather(local_result, root=0)

    if rank == 0:
        nonempty = [arr for arr in gathered if isinstance(arr, np.ndarray) and arr.size > 0]
        if not nonempty:
            raise RuntimeError(f"No results gathered for case: {case_name}")

        combined_results = np.concatenate(nonempty, axis=0)
        np.savetxt(output_file, combined_results, **np_savetxt_args)
        print(f"Wrote: {output_file}", flush=True)

    comm.Barrier()


def load_config(config_file: str) -> list[dict]:
    """
    Load YAML config file.

    Expected format:
    runs:
      - name: ...
        run_dir: ...
        var_name: ...
        output_file: ...
        bpfile: ...
        start_stack: 2
        end_stack: 366
    """
    with open(config_file, 'r') as f:
        cfg = yaml.safe_load(f)

    if isinstance(cfg, dict) and 'runs' in cfg:
        runs = cfg['runs']
    elif isinstance(cfg, list):
        runs = cfg
    else:
        raise ValueError("YAML config must be either a list of runs or a dict with key 'runs'.")

    if not isinstance(runs, list):
        raise ValueError("'runs' must be a list.")

    required = {'run_dir', 'var_name', 'output_file', 'bpfile', 'start_stack', 'end_stack'}
    for i, case in enumerate(runs):
        if not isinstance(case, dict):
            raise ValueError(f"Run #{i} is not a dictionary.")
        missing = required - set(case.keys())
        if missing:
            raise ValueError(f"Run #{i} is missing required keys: {sorted(missing)}")

    return runs


def run_from_config(config_file: str):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        print(f"Loading YAML config: {config_file}", flush=True)

    runs = load_config(config_file)

    for i, case in enumerate(runs):
        if rank == 0:
            print(f"\n--- Run {i + 1}/{len(runs)} ---", flush=True)
        process_one_case(case)

    if rank == 0:
        print("\nAll YAML runs finished.", flush=True)


def build_single_case_from_top_inputs() -> dict:
    return {
        'name': 'single_case_from_top_inputs',
        'run_dir': RUN_DIR,
        'var_name': VAR_NAME,
        'output_file': OUTPUT_FILE,
        'bpfile': BPFILE,
        'start_stack': start_stack,
        'end_stack': end_stack,
    }


def main():
    rank = MPI.COMM_WORLD.Get_rank()

    if CONFIG_FILE is None:
        if rank == 0:
            print("Running single case from top-of-script inputs.", flush=True)
        process_one_case(build_single_case_from_top_inputs())
    else:
        run_from_config(CONFIG_FILE)

    if rank == 0:
        print("done!")


if __name__ == '__main__':
    main()