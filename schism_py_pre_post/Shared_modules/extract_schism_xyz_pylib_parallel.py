import os
os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg", force=True)

import numpy as np
import xarray as xr
import dask
from dask import delayed
from dask.distributed import Client
from pylib import read_schism_output, save_schism_grid, read


# --------------- inputs ----------------
RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R35/'
output_files = [
    '/sciclone/schism10/feiye/STOFS3D-v8/O35/ST.transect.mississippi',
]
bpfile_list = [
    '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/mississippi.bp'
]

start_stack = 1
end_stack = 142
num_jobs = 20  # Adjust as needed

# ---------------------------------------
# format of the output file
np_savetxt_args = {'fmt': '%.4f', 'delimiter': ' ', 'newline': '\n'}


def process_stack_range(stack_range, bpfile):
    '''
    Process a range of stacks
    '''
    print(f"Processing stack range: {stack_range[0]} to {stack_range[-1]}\n")
    data = read_schism_output(run=RUN_DIR, varname=['elevation'], xyz=bpfile, stacks=stack_range)
    return np.c_[data.time, data.elevation.T]


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


def extract_transect(stack_range=None, bpfile=None):
    '''
    Test serial processing
    varname: string, e.g., 'salinity'
    stack_range: range of stacks to process, e.g., np.arange(1, 10)
    '''
    varname = ['salinity', 'temperature', 'zCoordinates']
    data = read_schism_output(
        run=RUN_DIR, varname=varname, xyz=bpfile, stacks=stack_range,
        fmt=1,  # transect xy and all levels
    )
    return data


def combine_results(results):
    '''
    Combine results from multiple Dask tasks into a single xarray Dataset
    '''
    variable_names = [info.split(':')[0].strip() for info in results[0].INFO]
    ds = xr.Dataset()
    for var_name in variable_names:
        time_dim = 0 if var_name == 'time' else 1
        var_combined = np.concatenate([eval(f"result.{var_name}") for result in results], axis=time_dim)

        if var_name == 'time':
            ds[var_name] = (('time',), var_combined)
        elif var_name in ['elevation', 'elev']:
            data = np.transpose(var_combined)  # transpose to (nt, npoints)
            ds[var_name] = (('time', 'npoints'), data)
        elif var_name in ['salinity', 'temperature', 'zCoordinates']:
            data = np.transpose(var_combined, (1, 0, 2))  # transpose to (nt, npoints, nz)
            ds[var_name] = (('time', 'npoints', 'nz'), data)
        else:
            raise ValueError(f"Unexpected number of dimensions: {var_combined.ndim} for variable {var_name}")

    return ds


def parallel_extract():
    '''
    Parallel extraction using Dask
    '''
    # Initialize the Dask client (This will allow you to monitor the task status in real-time)
    client = Client()

    # Split the work into chunks for Dask to process in parallel
    stack_ranges = np.array_split(np.arange(start_stack, end_stack + 1), num_jobs)
    # Use Dask's delayed to lazily compute each task (will not execute until explicitly computed)
    # tasks = [delayed(extract_transect)(sr, bpfile_list[0]) for sr in stack_ranges]
    tasks = [delayed(extract_transect)(sr, bpfile_list[0]) for sr in stack_ranges]
    client.close()

    # Run the tasks in parallel using Dask
    results = dask.compute(*tasks)
    results_ds = combine_results(results)
    # attach lon/lat from bpfile
    bp = read(bpfile_list[0])
    results_ds['lon'] = (('npoints',), bp.x)
    results_ds['lat'] = (('npoints',), bp.y)

    # write the results to file
    results_ds.to_netcdf(output_files[0] + '.nc')


def serial_extract():
    '''
    Test serial processing
    '''
    # attach lon/lat from bpfile
    bp = read(bpfile_list[0])

    results = extract_transect(
        stack_range=np.arange(start_stack, end_stack + 1),
        bpfile=bpfile_list[0]
    )
    results_ds = combine_results([results])
    results_ds['lon'] = (('npoints',), bp.x)
    results_ds['lat'] = (('npoints',), bp.y)

    # write the results to file
    results_ds.to_netcdf(output_files[0] + '.nc')


if __name__ == '__main__':
    # serial_extract()
    parallel_extract()
    print("done!")
