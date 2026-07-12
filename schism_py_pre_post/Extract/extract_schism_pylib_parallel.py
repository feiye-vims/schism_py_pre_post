#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import tempfile
from typing import Optional

from mpi4py import MPI
import numpy as np
import yaml

from pylib import read_schism_output


# --------------------------- default inputs ---------------------------
# Command-line arguments can select a different config or single-case mode.
CONFIG_FILE = '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/extract.yaml'

RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v7.3/F2018/'
VAR_NAME = 'elev'
OUTPUT_FILE = '/sciclone/schism10/feiye/STOFS3D-v7.3/F2018/elevation.stofs3d_atl_202501.dat'
BPFILE = '/sciclone/schism10/feiye/STOFS3D-v7.3/BPfiles/stofs3d_atl_202501.bp'
START_STACK = 2
END_STACK = 366

DEFAULT_OUTPUT_FORMAT = '%.4f'
NP_SAVETXT_ARGS = {'delimiter': ' ', 'newline': '\n'}
# ---------------------------------------------------------------------


def process_stack_range(
    run_dir: str,
    var_name: str,
    stack_range: np.ndarray,
    bpfile: str,
) -> Optional[np.ndarray]:
    """Extract one rank's stack range and return time-first tabular data."""
    if stack_range.size == 0:
        return None

    rank = MPI.COMM_WORLD.Get_rank()
    print(
        f"[rank {rank}] processing stacks {int(stack_range[0])}.."
        f"{int(stack_range[-1])} for run={run_dir}",
        flush=True,
    )

    data = read_schism_output(
        run=run_dir,
        varname=[var_name],
        xyz=bpfile,
        stacks=stack_range,
    )

    time = np.asarray(data.time)
    values = np.asarray(getattr(data, var_name))
    if time.ndim != 1:
        raise ValueError(f"Expected one-dimensional time data, got shape {time.shape}.")

    if values.ndim == 1 and values.size == time.size:
        values = values[:, np.newaxis]
    elif values.ndim == 2 and values.shape[1] == time.size:
        # pylib scalar output normally uses (point, time).
        values = values.T
    elif values.ndim == 2 and values.shape[0] == time.size:
        values = values
    else:
        raise ValueError(
            f"Variable {var_name!r} has shape {values.shape}; expected scalar data "
            f"with one time dimension of length {time.size}."
        )

    return np.column_stack((time, values))


def split_stacks(start_s: int, end_s: int, size: int, rank: int) -> np.ndarray:
    """Split the inclusive stack interval into contiguous rank-local chunks."""
    return np.array_split(np.arange(start_s, end_s + 1, dtype=int), size)[rank]


def validate_case(case: dict) -> dict:
    """Validate and normalize one case before any ranks start extraction."""
    required = {'run_dir', 'var_name', 'output_file', 'bpfile', 'start_stack', 'end_stack'}
    missing = required - set(case)
    if missing:
        raise ValueError(f"Case is missing required keys: {sorted(missing)}")

    normalized = dict(case)
    normalized['run_dir'] = str(Path(case['run_dir']).expanduser())
    normalized['bpfile'] = str(Path(case['bpfile']).expanduser())
    normalized['output_file'] = str(Path(case['output_file']).expanduser())
    normalized['var_name'] = str(case['var_name']).strip()
    normalized['start_stack'] = int(case['start_stack'])
    normalized['end_stack'] = int(case['end_stack'])

    if not Path(normalized['run_dir']).is_dir():
        raise FileNotFoundError(f"Run directory does not exist: {normalized['run_dir']}")
    if not Path(normalized['bpfile']).is_file():
        raise FileNotFoundError(f"BP file does not exist: {normalized['bpfile']}")
    if not normalized['var_name']:
        raise ValueError("var_name must not be empty.")
    if normalized['start_stack'] > normalized['end_stack']:
        raise ValueError(
            f"start_stack ({normalized['start_stack']}) exceeds "
            f"end_stack ({normalized['end_stack']})."
        )
    if Path(normalized['output_file']).is_dir():
        raise IsADirectoryError(f"Output path is a directory: {normalized['output_file']}")

    return normalized


def save_atomic(output_file: str, data: np.ndarray, fmt: str, overwrite: bool) -> None:
    """Write in the output directory and atomically publish the completed file."""
    destination = Path(output_file)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {destination}")

    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode='w',
            dir=destination.parent,
            prefix=f'.{destination.name}.',
            suffix='.tmp',
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            np.savetxt(temporary, data, fmt=fmt, **NP_SAVETXT_ARGS)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def process_one_case(
    case: dict,
    fmt_override: Optional[str] = None,
    overwrite_default: bool = True,
) -> None:
    """Extract, gather, and write one case with collective error handling."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    validation_error = None
    normalized_case = None
    if rank == 0:
        try:
            normalized_case = validate_case(case)
        except Exception as exc:  # Broadcast input failures instead of stranding ranks.
            validation_error = f"{type(exc).__name__}: {exc}"
    normalized_case, validation_error = comm.bcast(
        (normalized_case, validation_error), root=0
    )
    if validation_error is not None:
        raise RuntimeError(f"Case validation failed: {validation_error}")
    case = normalized_case

    run_dir = case['run_dir']
    var_name = case['var_name']
    output_file = case['output_file']
    bpfile = case['bpfile']
    start_s = case['start_stack']
    end_s = case['end_stack']
    case_name = case.get('name', 'unnamed_case')
    output_fmt = fmt_override or case.get('fmt', DEFAULT_OUTPUT_FORMAT)
    overwrite = overwrite_default

    if not isinstance(output_fmt, str) or not output_fmt:
        raise ValueError("Output format must be a nonempty NumPy savetxt format string.")

    if rank == 0:
        print(f"\n=== Processing case: {case_name} ===", flush=True)
        print(f"RUN_DIR    : {run_dir}", flush=True)
        print(f"VAR_NAME   : {var_name}", flush=True)
        print(f"BPFILE     : {bpfile}", flush=True)
        print(f"STACKS     : {start_s}..{end_s}", flush=True)
        print(f"OUTPUT     : {output_file}", flush=True)

    local_result = None
    local_error = None
    try:
        local_stacks = split_stacks(start_s, end_s, size, rank)
        local_result = process_stack_range(run_dir, var_name, local_stacks, bpfile)
    except Exception as exc:
        local_error = f"rank {rank}: {type(exc).__name__}: {exc}"

    extraction_errors = [error for error in comm.allgather(local_error) if error]
    if extraction_errors:
        raise RuntimeError("Extraction failed:\n" + '\n'.join(extraction_errors))

    gathered = comm.gather(local_result, root=0)
    write_error = None
    if rank == 0:
        try:
            nonempty = [array for array in gathered if array is not None]
            if not nonempty:
                raise RuntimeError(f"No results gathered for case: {case_name}")
            combined_results = np.concatenate(nonempty, axis=0)
            save_atomic(output_file, combined_results, output_fmt, overwrite)
            print(f"Wrote: {output_file}", flush=True)
        except Exception as exc:
            write_error = f"{type(exc).__name__}: {exc}"

    write_error = comm.bcast(write_error, root=0)
    if write_error is not None:
        raise RuntimeError(f"Output failed: {write_error}")


def load_config(config_file: str) -> list[dict]:
    """Load a YAML list of cases, optionally nested under a ``runs`` key."""
    with open(config_file, 'r', encoding='utf-8') as stream:
        cfg = yaml.safe_load(stream)

    if isinstance(cfg, dict) and 'runs' in cfg:
        runs = cfg['runs']
    elif isinstance(cfg, list):
        runs = cfg
    else:
        raise ValueError("YAML config must be a list of runs or a dict with key 'runs'.")

    if not isinstance(runs, list) or not runs:
        raise ValueError("'runs' must be a nonempty list.")
    if any(not isinstance(case, dict) for case in runs):
        raise ValueError("Every run must be a dictionary.")
    return runs


def select_runs(runs: list[dict], requested_names: Optional[list[str]]) -> list[dict]:
    """Select named runs while rejecting ambiguous names in the configuration."""
    runs_by_name = {}
    duplicate_names = set()
    for case in runs:
        name = case.get('name')
        if name is None:
            continue
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Every specified run name must be a nonempty string.")
        name = name.strip()
        if name in runs_by_name:
            duplicate_names.add(name)
        runs_by_name[name] = case

    if duplicate_names:
        raise ValueError(
            f"Duplicate run names in configuration: {sorted(duplicate_names)}"
        )
    if not requested_names:
        return runs

    unique_requested_names = list(dict.fromkeys(requested_names))
    missing_names = [name for name in unique_requested_names if name not in runs_by_name]
    if missing_names:
        available_names = sorted(runs_by_name)
        raise ValueError(
            f"Requested run names not found: {missing_names}. "
            f"Available names: {available_names}"
        )
    return [runs_by_name[name] for name in unique_requested_names]


def run_from_config(
    config_file: str,
    fmt_override: Optional[str],
    overwrite_default: bool,
    requested_names: Optional[list[str]] = None,
) -> None:
    """Load the configuration once on rank 0, then broadcast it to all ranks."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    runs = None
    load_error = None

    if rank == 0:
        print(f"Loading YAML config: {config_file}", flush=True)
        try:
            runs = select_runs(load_config(config_file), requested_names)
        except Exception as exc:
            load_error = f"{type(exc).__name__}: {exc}"

    runs, load_error = comm.bcast((runs, load_error), root=0)
    if load_error is not None:
        raise RuntimeError(f"Unable to load configuration: {load_error}")

    for index, case in enumerate(runs):
        if rank == 0:
            print(f"\n--- Run {index + 1}/{len(runs)} ---", flush=True)
        process_one_case(case, fmt_override, overwrite_default)

    if rank == 0:
        print("\nAll YAML runs finished.", flush=True)


def build_single_case_from_defaults() -> dict:
    return {
        'name': 'single_case_from_defaults',
        'run_dir': RUN_DIR,
        'var_name': VAR_NAME,
        'output_file': OUTPUT_FILE,
        'bpfile': BPFILE,
        'start_stack': START_STACK,
        'end_stack': END_STACK,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--config', help='YAML configuration file')
    mode.add_argument(
        '--single-case',
        action='store_true',
        help='use the single-case defaults in this script',
    )
    parser.add_argument(
        '--fmt',
        help=f"override the NumPy output format (default: {DEFAULT_OUTPUT_FORMAT})",
    )
    parser.add_argument(
        '--no-overwrite',
        action='store_true',
        help='fail if an output file already exists',
    )
    parser.add_argument(
        '--run',
        dest='run_names',
        action='append',
        metavar='NAME',
        help='process only this named YAML run; may be specified more than once',
    )
    args = parser.parse_args()
    if args.run_names and (args.single_case or (not args.config and CONFIG_FILE is None)):
        parser.error('--run requires YAML configuration (use --config if needed)')
    return args


def main() -> None:
    rank = MPI.COMM_WORLD.Get_rank()
    args = parse_args()
    config_file = None if args.single_case else (args.config or CONFIG_FILE)

    if config_file is None:
        if rank == 0:
            print("Running single case from script defaults.", flush=True)
        process_one_case(
            build_single_case_from_defaults(),
            fmt_override=args.fmt,
            overwrite_default=not args.no_overwrite,
        )
    else:
        run_from_config(
            config_file,
            fmt_override=args.fmt,
            overwrite_default=not args.no_overwrite,
            requested_names=args.run_names,
        )

    if rank == 0:
        print("done!", flush=True)


if __name__ == '__main__':
    main()
