import xarray as xr
from pathlib import Path
import numpy as np
from pylib_essentials.schism_file import cread_schism_hgrid

# dictionary of variable names and their residing schism output files
var_dict = {
    'elevation':
        {"file_prefix": "out2d"},
}

def max_var(output_dir:Path, var_str:str, stacks=[]):
    """Return the maximum value of a variable
    from specified stacks of schism outputs."""

    if stacks == []:  # if no stacks are specified, use all files in the output directory
        files = [f for f in output_dir.glob(f'*{var_dict[var_str]["file_prefix"]}*.nc')]
    else:
        files = [f'{output_dir}/{var_dict[var_str]["file_prefix"]}_{stack}.nc' for stack in stacks]

    ds = xr.open_mfdataset(files, combine='by_coords')
    var = ds[var_str]
    var = var.max(dim='time')
    return var.values

if __name__ == "__main__":
    # example usage
    output_dir = Path('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R14b/outputs/')
    stacks = np.arange(31, 32)
    max_elev = max_var(output_dir, "elevation", stacks)

     # output as *.2dm
    hgrid_obj = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R14b/hgrid.gr3')
    original_dp = hgrid_obj.dp.copy()

    hgrid_obj.dp = max_elev  # max elevation
    hgrid_obj.dp[original_dp + max_elev < 0.01] = -9999 # mask out shallow inundation
    hgrid_obj.grd2sms(f'{output_dir}/max_elev_stack31.2dm')

    hgrid_obj.dp = max_elev + hgrid_obj.dp
    hgrid_obj.dp[original_dp + max_elev < 0.01] = -9999 # mask out shallow inundation
    hgrid_obj.grd2sms(f'{output_dir}/max_inun_stack31.2dm')

    pass
