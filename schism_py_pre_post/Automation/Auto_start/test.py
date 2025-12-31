#!/usr/bin/env python3

import os
import time
import subprocess
import glob
import re
import time
from datetime import datetime

#----------input---------------------------------
#(1) Copy this script, run_test and run_comb into rundir

#(2) Inside rundir: prep inputs as before
#    The "hotstart.nc" (if it exists under ihot=1 or 2) will be overwritten by this script by symlinks ("ln -sf"); 
#    if you want to keep it, "mv hotstart.nc hotstart.nc.0" then "ln -s hotstart.nc.0 hotstart.nc"

#(3) Under the run dir,
#      "python auto_hotstart.py >& scrn.out"
#    or
#      "./auto_hotstart.py > scrn.out &"
#      "tail -f scrn.out"
#    The latter prints jobs status and run speeds to the screen continuously.

#    This can be done at any stage of the run, for example:
#    * Before a run is launched, i.e., the script will submit the initial job and monitor its progress;
#    * After the run is launched manually;
#    * After the run is interrupted, e.g., due to time limit or manual scancel

#    Note: the script will use the last part of the current dir as runid, e.g., the runid of "RUN13a" or "R13a" will be "13a"
#          , make sure you don't have duplicate run ids.

#(4) Sometimes a run can hang, e.g., on Hercules, the script will automatically detect this and interrupt the run then relaunch from
#    the nearest hotstarting point

job_scheduler = 'slurm' # 'pbs' or 'slurm'

rundir = os.getcwd()  # use os.getcwd if launch from rundir; otherwise specify a string yourself
last_stack = None  # if None, the script will try to find the last stack number in the file "param.nml"
                   # make sure the run can finish the specified rnday in param.nml (i.e., the forcing covers the whole period);
                   # otherwise, change the "rnday" in param.nml or specify another number here

#----------end input---------------------------------


# ---------------------  embedded functions -----------------------

def Replace_string_in_file(fname, str_orig_pattern, str_replace):
    if '~' in fname:
        fname = fname.replace('~', os.path.expanduser('~'))
    fname = os.path.abspath(fname)
    with open(fname, "rt") as fin:
        with open("tmp.txt", "wt") as fout:
            for line in fin:
                fout.write(re.sub(str_orig_pattern, str_replace, line))
    os.system(f"mv tmp.txt {fname}")


def ReplaceJobName(fname, job_name, job_scheduler='slurm'):
    import fileinput

    if job_scheduler == 'slurm':
        pattern = r"(SBATCH\s+-J\s+)(\S+)"
    else:
        raise Exception('job_scheduler must be either "slurm"; pbs not implemented yet')

    replacement = rf'SBATCH -J {job_name}'

    # Use fileinput to edit the file in place
    match_found = False
    with fileinput.FileInput(fname, inplace=True) as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                match_found = True
            modified_line = re.sub(pattern, replacement, line)
            print(modified_line, end='')
    
    if not match_found:
        raise Exception(f'Job name specification not found in {fname}')

def Get_hotstart_step(run_out_dir):
    if '~' in run_out_dir:
        run_out_dir = run_out_dir.replace('~', os.path.expanduser('~'))
    run_out_dir = os.path.abspath(run_out_dir)

    hot_files = glob.glob(f"{run_out_dir}/hotstart_000000_*.nc")

    hot_steps = []
    for hot_file in hot_files:
        sub_str = hot_file.split("_")
        sub_str = sub_str[-1].split(".")
        hot_steps.append(int(sub_str[0]))

    hot_steps.sort()
    return hot_steps

def Get_var_from_file(fname, var_name, reverse_search=False):
    if '~' in fname:
        fname = fname.replace('~', os.path.expanduser('~'))
    fname = os.path.abspath(fname)


    pattern = fr'{var_name}\s*=\s*(\d+)'
    with open(fname, "rt") as fin:
        lines = fin.readlines()  # Read all lines into a list

    if reverse_search:
        lines =[line for line in reversed(lines)]

    for line in lines:
        match = re.search(pattern, line)
        if match:
            var_value = match.group(1)  # Extract the matched group
            return var_value

    print(f'Variable {var_name} not found in {fname}')
    return None

# --------------------- end embedded functions -----------------------

rundir = './'

Replace_string_in_file(f'{rundir}/run_comb', '-i 0000', f'-i 12345')


Replace_string_in_file(f'{rundir}/param.nml', r'ihot\s*=\s*\d+', 'ihot = 2')
print(f'Task completed.', flush=True)
