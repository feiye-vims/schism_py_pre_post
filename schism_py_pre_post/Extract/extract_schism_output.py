"""Helpers for extracting SCHISM station time series on demand.

This module validates raw SCHISM stack layouts, builds the YAML/JSON config
consumed by ``extract_schism_pylib_parallel.py``, prompts before running MPI,
and rechecks the expected model file before extraction starts.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_EXTRACT_SCRIPT = "Extract.extract_schism_pylib_parallel"
TRAILING_STACK_SIZE_TOLERANCE = 0.01
TRAILING_STACK_SIZE_SIGMA = 3.0


def stack_files(raw_outputs_dir: Path, stack_prefix: str) -> dict[int, Path]:
    pattern = re.compile(rf"^{re.escape(stack_prefix)}_(\d+)\.nc$")
    stacks = {}
    for path in raw_outputs_dir.iterdir():
        match = pattern.match(path.name)
        if match:
            stacks[int(match.group(1))] = path
    return dict(sorted(stacks.items()))


def format_stack_list(stacks: list[int]) -> str:
    if not stacks:
        return "none"
    ranges = []
    start = previous = stacks[0]
    for stack in stacks[1:]:
        if stack == previous + 1:
            previous = stack
            continue
        ranges.append(f"{start}" if start == previous else f"{start}-{previous}")
        start = previous = stack
    ranges.append(f"{start}" if start == previous else f"{start}-{previous}")
    return ", ".join(ranges)


def stack_prefix_candidates(var_name: str | None) -> list[str]:
    if not var_name:
        return ["schout", "out2d"]
    normalized = str(var_name).strip().lower()
    aliases = {
        "elev": ["out2d", "schout"],
        "elevation": ["out2d", "schout"],
        "depth": ["out2d", "schout"],
        "salinity": ["salinity"],
        "temperature": ["temperature"],
        "horizontalvelx": ["horizontalVelX"],
        "horizontalvely": ["horizontalVelY"],
        "zcoordinates": ["zCoordinates"],
    }
    return aliases.get(normalized, [var_name])


def infer_stack_range(
    output_dir: str | os.PathLike[str],
    var_name: str | None = None,
    trailing_stack_size_tolerance: float = TRAILING_STACK_SIZE_TOLERANCE,
    trailing_stack_size_sigma: float = TRAILING_STACK_SIZE_SIGMA,
) -> dict[str, Any]:
    output_dir = Path(output_dir)
    raw_outputs_dir = output_dir / "outputs"
    if not output_dir.is_dir():
        raise FileNotFoundError(f"Run output_dir does not exist: {output_dir}")
    if not raw_outputs_dir.is_dir():
        raise FileNotFoundError(f"Expected raw SCHISM outputs directory does not exist: {raw_outputs_dir}")

    stack_sets = {
        prefix: stack_files(raw_outputs_dir, prefix)
        for prefix in stack_prefix_candidates(var_name)
    }
    present = {name: files for name, files in stack_sets.items() if files}
    if not present:
        expected_patterns = " or ".join(f"{prefix}_*.nc" for prefix in stack_sets)
        raise FileNotFoundError(
            f"No SCHISM stack files found in {raw_outputs_dir}; expected {expected_patterns}"
        )
    if len(present) != 1:
        found_patterns = ", ".join(f"{prefix}_*.nc" for prefix in present)
        raise ValueError(
            f"Ambiguous SCHISM stack layout in {raw_outputs_dir}; found {found_patterns}"
        )

    stack_type, stack_file_map = next(iter(present.items()))
    stacks = list(stack_file_map)
    expected = list(range(stacks[0], stacks[-1] + 1))
    if stacks != expected:
        missing = sorted(set(expected) - set(stacks))
        preview = missing[:10]
        suffix = "..." if len(missing) > len(preview) else ""
        raise ValueError(
            f"Non-consecutive {stack_type}_*.nc stacks in {raw_outputs_dir}; missing {preview}{suffix}"
        )

    excluded_stacks = []
    included_stacks = list(stacks)
    stack_sizes = {stack: stack_file_map[stack].stat().st_size for stack in stacks}
    if len(stacks) > 1:
        last_stack = stacks[-1]
        reference_sizes = [stack_sizes[stack] for stack in stacks[:-1]]
        reference_mean = sum(reference_sizes) / len(reference_sizes)
        if len(reference_sizes) > 1:
            reference_std = (
                sum((size - reference_mean) ** 2 for size in reference_sizes)
                / (len(reference_sizes) - 1)
            ) ** 0.5
        else:
            reference_std = 0.0
        lower_bound = reference_mean * (1.0 - trailing_stack_size_tolerance)
        if reference_std > 0.0:
            lower_bound = reference_mean - trailing_stack_size_sigma * reference_std
        if stack_sizes[last_stack] < lower_bound:
            excluded_stacks.append(last_stack)
            included_stacks = stacks[:-1]

    if not included_stacks:
        raise ValueError(f"No complete {stack_type}_*.nc stacks remain after size checks in {raw_outputs_dir}")

    expected_included = list(range(included_stacks[0], included_stacks[-1] + 1))
    if included_stacks != expected_included:
        raise ValueError(
            f"Included {stack_type}_*.nc stacks are not consecutive after exclusions in {raw_outputs_dir}: "
            f"{format_stack_list(included_stacks)}"
        )

    return {
        "run_dir": output_dir,
        "raw_outputs_dir": raw_outputs_dir,
        "start_stack": included_stacks[0],
        "end_stack": included_stacks[-1],
        "stack_type": stack_type,
        "included_stacks": included_stacks,
        "excluded_stacks": excluded_stacks,
        "stack_sizes": stack_sizes,
        "size_reference_mean": reference_mean if len(stacks) > 1 else stack_sizes[stacks[0]],
        "size_reference_std": reference_std if len(stacks) > 1 else 0.0,
        "size_lower_bound": lower_bound if len(stacks) > 1 else stack_sizes[stacks[0]],
    }


def extraction_case_from_run_spec(
    run_spec: dict[str, Any],
    trailing_stack_size_tolerance: float = TRAILING_STACK_SIZE_TOLERANCE,
) -> dict[str, Any]:
    extract = run_spec.get("extract")
    if extract is None:
        extract = {}
    if not isinstance(extract, dict):
        raise TypeError(f"extract block for {run_spec['name']} must be a dictionary")
    obsolete_keys = {"run_dir", "output_file", "start_stack", "end_stack"} & set(extract)
    if obsolete_keys:
        raise ValueError(
            f"extract block for {run_spec['name']} must not define {sorted(obsolete_keys)}; "
            "these are inferred from output_dir and elev_out_file"
        )

    if "output_dir" not in run_spec:
        raise FileNotFoundError(
            f"Model file does not exist and output_dir is not configured for {run_spec['name']}: "
            f"{run_spec['elev_out_file']}"
        )

    case = dict(extract)
    case.setdefault("name", run_spec["name"])
    case.setdefault("var_name", "elevation")
    stack_info = infer_stack_range(
        run_spec["output_dir"],
        var_name=case["var_name"],
        trailing_stack_size_tolerance=trailing_stack_size_tolerance,
    )
    case.setdefault("run_dir", str(stack_info["run_dir"]))
    case.setdefault("output_file", run_spec["elev_out_file"])
    case.setdefault("bpfile", run_spec["station_bp_file"])
    case.setdefault("start_stack", stack_info["start_stack"])
    case.setdefault("end_stack", stack_info["end_stack"])
    case["stack_type"] = stack_info["stack_type"]
    case["_stack_info"] = stack_info

    required = {"run_dir", "var_name", "output_file", "bpfile", "start_stack", "end_stack"}
    missing = required - set(case)
    if missing:
        raise ValueError(f"Extract block for {run_spec['name']} is missing required keys: {sorted(missing)}")

    return case


def write_extract_config(case: dict[str, Any], config_dir: str | os.PathLike[str]) -> Path:
    config_dir = Path(config_dir)
    config_dir.mkdir(parents=True, exist_ok=True)
    config_file = config_dir / f"extract_{case['name']}.yaml"
    with open(config_file, "w", encoding="utf-8") as f:
        json.dump({"runs": [case]}, f, indent=2)
        f.write("\n")
    return config_file


def model_file_ready(elev_out_file: str | os.PathLike[str]) -> bool:
    elev_out_file = Path(elev_out_file)
    return elev_out_file.is_file() and elev_out_file.stat().st_size > 0


def confirm_extraction(
    run_name: str,
    command: list[str],
    elev_out_file: Path,
    config_file: Path,
    stack_info: dict[str, Any],
) -> None:
    print(f"\nModel file not found for {run_name}: {elev_out_file}")
    print(f"Extraction config: {config_file}")
    print(f"Raw outputs: {stack_info['raw_outputs_dir']}")
    print(f"Stack type: {stack_info['stack_type']}")
    print(f"Stacks to combine: {format_stack_list(stack_info['included_stacks'])}")
    print(f"Stacks excluded: {format_stack_list(stack_info['excluded_stacks'])}")
    if stack_info["excluded_stacks"]:
        last_stack = stack_info["excluded_stacks"][-1]
        previous_stack = stack_info["included_stacks"][-1]
        sizes = stack_info["stack_sizes"]
        print(
            f"Excluded trailing stack {last_stack} because its file size "
            f"({sizes[last_stack]} bytes) is below mean - {TRAILING_STACK_SIZE_SIGMA:g} std "
            f"({stack_info['size_lower_bound']:.0f} bytes). Previous included stack "
            f"{previous_stack} is {sizes[previous_stack]} bytes."
        )
    print("Extraction command:")
    print(shlex.join(command))
    answer = input("Proceed with extraction? [y/N]: ").strip().lower()
    if answer not in {"y", "yes"}:
        raise RuntimeError(f"Extraction cancelled for {run_name}")


def ensure_model_file(
    run_spec: dict[str, Any],
    extract_missing: bool = True,
    extract_mpi_np: int = 4,
    extract_script: str | os.PathLike[str] = DEFAULT_EXTRACT_SCRIPT,
    extract_config_dir: str | os.PathLike[str] = ".",
    prompt: bool = True,
) -> None:
    elev_out_file = Path(run_spec["elev_out_file"])
    if model_file_ready(elev_out_file):
        return

    if not extract_missing:
        raise FileNotFoundError(f"Model file does not exist: {elev_out_file}")

    case = extraction_case_from_run_spec(run_spec)
    stack_info = case.pop("_stack_info")
    case["output_file"] = str(elev_out_file)
    config_file = write_extract_config(case, extract_config_dir)
    nproc = int(case.pop("nproc", extract_mpi_np))
    fmt = case.pop("fmt", None)

    extract_target = str(extract_script)
    python_args = (
        ["python", extract_target]
        if extract_target.endswith(".py") or "/" in extract_target
        else ["python", "-m", extract_target]
    )
    command = [
        "mpiexec", "-n", str(nproc), *python_args,
        "--config", str(config_file), "--run", case["name"],
    ]
    if fmt:
        command.extend(["--fmt", str(fmt)])

    if prompt:
        confirm_extraction(run_spec["name"], command, elev_out_file, config_file, stack_info)
    else:
        print("Extraction command:")
        print(shlex.join(command))

    if model_file_ready(elev_out_file):
        print(f"Model file appeared before extraction started; using existing file: {elev_out_file}")
        return

    subprocess.run(command, check=True)

    if not model_file_ready(elev_out_file):
        raise RuntimeError(f"Extraction finished but output is missing or empty: {elev_out_file}")


def _load_run_spec(path: str | os.PathLike[str]) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_spec_json", help="JSON file containing one run spec")
    parser.add_argument("--no-extract-missing", action="store_true")
    parser.add_argument("--extract-mpi-np", type=int, default=4)
    parser.add_argument("--extract-script", default=str(DEFAULT_EXTRACT_SCRIPT))
    parser.add_argument("--extract-config-dir", default=".")
    parser.add_argument("--yes", action="store_true", help="run extraction without interactive confirmation")
    args = parser.parse_args()

    ensure_model_file(
        _load_run_spec(args.run_spec_json),
        extract_missing=not args.no_extract_missing,
        extract_mpi_np=args.extract_mpi_np,
        extract_script=args.extract_script,
        extract_config_dir=args.extract_config_dir,
        prompt=not args.yes,
    )


if __name__ == "__main__":
    main()
