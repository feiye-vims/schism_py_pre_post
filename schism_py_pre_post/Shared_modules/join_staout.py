"""
For EnOI runs, join the station output files from individual runs
and link all schout_*.nc files into a single directory.
"""
import os
import numpy as np
from glob import glob
import re

rundir = '/sciclone/schism10/feiye/STOFS3D-v7.3/r2009/'
days_per_stack = 1

sub_runs = sorted(glob(f'{rundir}/start*/'))
aggregate_run = f'{rundir}/outputs/'
os.makedirs(aggregate_run, exist_ok=True)

# There should be a one-stack overlap between consecutive sub-runs,
# the last sub-run contains all data, but with later stack overwriting previous ones at overlaps
# , so we only need to replace the overlapping entries with the previous runs' data,
# , i.e., the last stack of each previous run.
data = np.loadtxt(f'{sub_runs[-1]}/schism_001/outputs/staout_1')
dt = data[1, 0] - data[0, 0]  # time step in seconds

previous_run_last_stack = 1
for i, sub_run in enumerate(sub_runs):
    sub_data = np.loadtxt(f'{sub_run}/schism_001/outputs/staout_1')
    stack_files = sorted([
        f for f in glob(f'{sub_run}/schism_*/outputs/schout_*.nc')
        if re.search(r'schout_\d+\.nc$', f)
    ], key=lambda x: int(re.search(r'schout_(\d+)\.nc$', x).group(1)))

    first_stack = int(re.search(r'schout_(\d+)\.nc$', stack_files[0]).group(1))
    if previous_run_last_stack != first_stack:
        raise ValueError(f"No overlap in stacks between runs: {previous_run_last_stack} vs {last_stack} under {sub_run}")

    last_stack = int(re.search(r'schout_(\d+)\.nc$', stack_files[-1]).group(1))
    last_stack_entries = np.arange(
        int((last_stack - 1) * days_per_stack * 86400 / dt),
        int(last_stack * days_per_stack * 86400 / dt)
    )

    # overwrite last stack in the main data
    data[last_stack_entries, 1:] = sub_data[last_stack_entries, 1:]

    # link schout files
    for f in stack_files:
        os.system(f'ln -s {f} {aggregate_run}/')  # don't use -sf to priroritize previous stack files

    previous_run_last_stack = last_stack


np.savetxt(f'{rundir}/staout_1', data, fmt='%.4f', delimiter=' ', newline='\n')
pass
