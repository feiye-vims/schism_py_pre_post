import os
import numpy as np
import dask
from dask import delayed
from dask.distributed import Client
from pylib import read_schism_output, save_schism_grid
from pylib_experimental.schism_file import TimeHistory

# --------------- inputs ----------------
RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R15a_v7/'
output_files = [
    '/sciclone/schism10/feiye/STOFS3D-v8/O15a_v7/outputs/elevation.coops_test.dat'
]
bpfile_list = [
    '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/coops_test.bp'
]

start_stack = 1
end_stack = 396
num_jobs = 12  # Adjust as needed

# ---------------------------------------
# format of the output file
np_savetxt_args = {'fmt': '%.4f', 'delimiter': ' ', 'newline': '\n'}


def process_stack_range(stack_range, bpfile):
    '''
    Process a range of stacks
    '''
    print(f"Processing stack range: {stack_range[0]} to {stack_range[-1]}\n")
    data = read_schism_output(run=RUN_DIR, varname=['elevation'], xyz=bpfile, stacks=stack_range)
    th = TimeHistory(data_array=np.c_[data.time, data.elevation.T])
    return np.c_[th.time + stack_range[0] - 1 + 1/24, th.data]


def test(num_range):
    '''
    Simple test function
    '''
    return num_range ** 2


def sample_save_grid():
    '''Sample function to save the grid'''

    rundir = '/sciclone/schism10/feiye/STOFS3D-v8/R13p_v7/'
    save_schism_grid(path=rundir, method=1)  # save full geometry)
    os.system(f"mv grid.npz {rundir}")


def main():
    '''
    Main function
    '''
    # Initialize the Dask client (This will allow you to monitor the task status in real-time)
    client = Client()

    # Split the work into chunks for Dask to process in parallel
    stack_ranges = np.array_split(np.arange(start_stack, end_stack + 1), num_jobs)

    # Use Dask's delayed to lazily compute each task (will not execute until explicitly computed)
    tasks = [delayed(process_stack_range)(sr, bpfile_list[0]) for sr in stack_ranges]

    # Run the tasks in parallel using Dask
    results = dask.compute(*tasks)
    combined_results = np.concatenate(results)

    # write the results to file
    np.savetxt(output_files[0], combined_results, **np_savetxt_args)

    client.close()
    return combined_results


if __name__ == '__main__':
    # sample_save_grid()
    main()
    print("done!")
