#!/usr/bin/env python

"""
Compare two netcdf files,
and print out the differences.
"""


import xarray as xr
import numpy as np
import argparse


# parse arguments
parser = argparse.ArgumentParser(description='Compare two netcdf files')
parser.add_argument('nc1', type=str, help='Path to the first netcdf file')
parser.add_argument('nc2', type=str, help='Path to the first netcdf file')
args = parser.parse_args()

nc1 = xr.open_dataset(args.nc1)
nc2 = xr.open_dataset(args.nc2)

# compare two nc files
# nc1 = xr.open_dataset(
#     '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/extract/schout_UV4.5m_1.nc')
# nc2 = xr.open_dataset(
#     '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/extract_copy/schout_UV4.5m_1.nc')

if sorted(list(nc1.variables.keys())) != sorted(list(nc2.variables.keys())):
    raise ValueError('Variables are not the same')
print('Variable lists are the same')

var_list = list(nc1.variables.keys())
print(f'Comparing variables values for {var_list}')

for key in var_list:
    dtype = nc1[key].dtype
    print(f'Comparing {key}, dtype: {dtype}')
    if np.issubdtype(dtype, np.floating):
        is_equal = np.allclose(np.array(nc1[key]), np.array(nc2[key]), equal_nan=True, atol=1e-6, rtol=1e-6)
    else:
        is_equal = nc1[key].equals(nc2[key])
    if not is_equal:
        raise ValueError(f'{key} is not equal')

print('All variables are the same')