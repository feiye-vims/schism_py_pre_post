import os
from glob import glob
import multiprocessing
from pylib import convert_dem_format

fnames = sorted(glob('/sciclone/schism10/Hgrid_projects/DEMs/CONED_2022/lonlat/*.tif'))

def tif2npz(arg_list):
    print(f"converting {arg_list[0]} to {arg_list[1]}")
    input_fname = arg_list[0]
    output_fname = arg_list[1]
    fmt = arg_list[2]
    convert_dem_format(input_fname, output_fname, fmt)

arg_list = [(fname, fname.replace('lonlat', 'npz').replace('tif', 'npz'), 1) for fname in fnames]
# tif2npz(arg_list[0])
with multiprocessing.Pool(processes=os.cpu_count()) as pool:
    pool.map(tif2npz, arg_list)

pass
