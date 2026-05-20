"""
For EnOI runs, join the yearly station output files into a single file
"""
import os
import numpy as np

proj_dir = '/sciclone/schism10/feiye/STOFS3D-v7.3/'
years = range(2000, 2026)
run_prefix = 'F'  # e.g., F2000, F2001, ..., F2025
aggregate_run = f'{proj_dir}/Total_Filtered/outputs/'
time_series_file = 'elevation.stofs3d_atl_202501.dat'

os.makedirs(aggregate_run, exist_ok=True)

total_base_time = np.datetime64(f'{years[0]}-01-01T00:00:00').astype('datetime64[s]').astype(float)
last_indices = None  # for diagnostic purpose
for i, year in enumerate(years):
    run = f'{proj_dir}/{run_prefix}{year}/'
    data = np.loadtxt(f'{run}/{time_series_file}')

    dts = sorted(np.diff(data[:, 0]))
    dt = np.mean(dts[int(len(dts)*0.4) : int(len(dts)*0.6)])  # use the middle 20% of time steps to avoid outliers
    if dt < 0.1:
        print(f"dt < 0.1, so probably in days, converting to seconds.")
        dt = int(dt * 24) * 3600  # round to hours, then convert to seconds
        data[:, 0] = np.round(data[:, 0] * 24) * 3600
        print(f"inferred dt: {dt} seconds.")

    base_time = np.datetime64(f'{year-1}-12-31T00:00:00').astype('datetime64[s]').astype(float)
    start_time = np.datetime64(f'{year}-01-01T00:00:00').astype('datetime64[s]').astype(float)
    end_time = np.datetime64(f'{year+1}-01-01T00:00:00').astype('datetime64[s]').astype(float)
    data_time = data[:, 0] + base_time

    if i == 0:  # first year, include start time
        mask = (data_time >= start_time) & (data_time <= end_time)
        # initialize total data array using first year's dt
        n_total_steps = int((np.datetime64(f'{years[-1]+1}-01-01T00:00:00').astype('datetime64[s]').astype(float) - total_base_time) / dt) + 1
        total_data = np.zeros((n_total_steps, data.shape[1]))
        total_data[:] = np.nan  # initialize with NaN
        print(f'Total data array initialized with {n_total_steps} time steps and {data.shape[1]-1} stations.')
    else:
        mask = (data_time > start_time) & (data_time <= end_time)

    total_data_indices = ((data_time[mask] - total_base_time) / dt).astype(int)

    if i > 0:
        if last_indices[-1] + 1 != total_data_indices[0]:
            raise ValueError(f"Data time steps are not continuous between years {year-1} and {year}.")

    total_data[total_data_indices, 0] = (data_time[mask] - total_base_time)
    total_data[total_data_indices, 1:] = data[mask, 1:]

    if i > 0:
        if total_data[total_data_indices[0]][0] != total_data[last_indices[-1]][0] + dt:
            raise ValueError(f"Data time steps are not continuous between years {year-1} and {year}.")

    last_indices = total_data_indices.copy()  # for diagnostic purpose

    print(f'Year: {year}, data points included: {np.sum(mask)}; inserted at indices: {total_data_indices[0]} to {total_data_indices[-1]}')

    pass

# remove any remaining NaN rows (if any)
nan_rows = np.isnan(total_data[:, 0])
total_data = total_data[~nan_rows, :]
if np.isnan(total_data).any():
    raise ValueError("There are still NaN values in the total data array after processing all years.")

# write the total data to file
np.savetxt(f'{aggregate_run}/{time_series_file}', total_data, fmt='%.4f', delimiter=' ', newline='\n')
