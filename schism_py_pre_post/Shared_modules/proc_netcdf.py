import os
import xarray as xr
import numpy as np
from copy import deepcopy
from pathlib import Path
from glob import glob


def fix_bottom_junk_nd(x, lo=-1e6, hi=1e6, require_all_components=True):
    """
    Replace invalid bottom layers (z=0 upward) with the first valid upper layer value.

    x: (nt, np, nz, ncomp) where z is bottom->surface
    lo/hi: validity bounds (adjust for your variable)
    require_all_components:
        - True: a layer is valid only if *all* components are valid (good for velocity u/v)
        - False: valid if *any* component is valid (less common)
    """
    assert x.ndim == 4, "Expected (nt, np, nz, ncomp)"
    nt, np_, nz, ncomp = x.shape

    # invalid per component
    invalid_c = (~np.isfinite(x)) | (x < lo) | (x > hi)      # (nt,np,nz,ncomp)
    valid_c = ~invalid_c

    if not invalid_c.any():
        print("No invalid values found.")
        return x

    invalid_np = invalid_c.any(axis=(0, 2, 3))   # shape: (np,)
    bad_nodes = np.where(invalid_np)[0]
    print("Nodes with invalid values:", bad_nodes)
    print(f"Number of bad nodes: {bad_nodes.size} out of {np_}")

    # reduce across components to decide if the *layer* is valid
    if require_all_components:
        valid_layer = valid_c.all(axis=3)                    # (nt,np,nz)
    else:
        valid_layer = valid_c.any(axis=3)                    # (nt,np,nz)

    has_valid = valid_layer.any(axis=2)                      # (nt,np)
    if not np.all(has_valid):
        raise ValueError("Some points have no valid layers at all.")

    first_valid_k = np.argmax(valid_layer, axis=2)           # (nt,np); 0 if all invalid

    # first valid values for all components at that k
    first_valid_val = np.take_along_axis(
        x, first_valid_k[..., None, None], axis=2
    )                                                        # (nt,np,1,ncomp)

    z = np.arange(nz)[None, None, :, None]                   # (1,1,nz,1)

    # Fill bottom junk layers (< first valid k).
    # You can choose whether to fill only layers that are invalid (recommended),
    # or overwrite everything below k. We'll fill only invalid entries.
    fill_mask = (z < first_valid_k[..., None, None]) & invalid_c  # (nt,np,nz,ncomp)

    x_fixed = x.copy()
    x_fixed[fill_mask] = np.broadcast_to(first_valid_val, x.shape)[fill_mask]

    return x_fixed

    
for ncfile in ['SAL_3D.th.nc', 'TEM_3D.th.nc', 'uv3D.th.nc', 'elev2D.th.nc']:
    # fill below-bottom junk values with bottom values
    if ncfile in ['SAL_3D.th.nc', 'TEM_3D.th.nc', 'uv3D.th.nc']:
        var = np.array(ds['time_series'][:])
        var_fixed = fix_bottom_junk_nd(var, lo=-1e6, hi=1e6, require_all_components=True)
        ds['time_series'] = (ds['time_series'].dims, var_fixed)
    elif ncfile == 'elev2D.th.nc':
        var = np.array(ds['time_series'][:])
        var += 0.2  # add an offset
        ds['time_series'] = (ds['time_series'].dims, var)

    ds.to_netcdf(f'{wdir}/{Path(ncfile).stem}_fix_bottom.nc')

    pass


# output_dir = '/sciclone/home/feiye/s1/STOFS3D-v7.3/R21g/outputs'
# ncfile_types = ['out2d', 'salinity', 'temperature', 'horizontalVelX', 'horizontalVelY', 'zCoordinates']
# for ncfile_type in ncfile_types:
#     ncfiles = sorted(glob(f'{output_dir}/{ncfile_type}_*.nc'))
#     for i, ncfile in enumerate(ncfiles):
#         os.symlink(ncfile, f'{wdir}/{ncfile_type}_{i+1}.nc')