import numpy as np
import re
import copy
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path
from pylib_experimental.schism_file import cread_schism_hgrid
from pylib import schism_grid


default_caxis = {
    'elevation': [-2, 2],
    'zeta': [-2, 2],
    'salinity': [0, 35],
    'temperature': [0, 35],
    'disturbance': [0, 2],
}


def match_out2d_fname(fname):
    """
    return corresponding out2d fname for a given netcdf file
    """
    pattern = r'(\w+)_(\d+)\.nc'  # \w+ captures any word before the number
    match = re.search(pattern, fname)
    
    if match:
        # Extract the prefix (group 1) and the number (group 2)
        prefix = match.group(1)
        number = match.group(2)
        
        out2d_fname = fname.replace(f'{prefix}_{number}.nc', f'out2d_{number}.nc')
    else:
        raise ValueError(f"Invalid filename {fname}")
    
    return out2d_fname
        

def read_schism_data(gd_fname, fnames, var_name, dry_mask_var, layer):
    """
    Read Schism grid and variable data.
    """
    print('Reading data ...')
    gd = cread_schism_hgrid(gd_fname)
    my_nc = xr.open_mfdataset(fnames)
    var = np.array(my_nc[var_name])
    my_nc.close()

    if dry_mask_var is not None or layer == "b":
        # read dry mask from out2d_*.nc instead of the original fnames
        fnames = [match_out2d_fname(fname) for fname in fnames]
        my_nc = xr.open_mfdataset(fnames)
        dry_mask = np.array(my_nc[dry_mask_var])
        kbp = np.array(my_nc['bottom_index_node'])
    else:
        dry_mask = np.zeros_like(var)
        kbp = None
    my_nc.close()

    return gd, var, dry_mask, kbp


def process_variable(var, dry_mask, it=-1, layer=-1, kbp=None):
    """
    Process the variable (e.g., mask dry nodes, extract surface or volume values).

    default input:
        it = -1: last time step
        layer = -1: surface value; "b" for bottom, "s" for surface, or layer number
        kbp = None: bottom layer index (needed for layer="b")
    """
    for i, _ in enumerate(var):
        var[i, dry_mask[i, :].astype(bool)] = np.nan

    if layer == "b":  # bottom value, needs kbp from vgrid
        if len(var.shape) == 3:
            value = var[it, :, :]
            value = value[np.arange(value.shape[0]), kbp]  # bottom index varies with nodes
        elif len(var.shape) == 2:
            value = var[it, :]
    elif layer == "s":  # surface value
        if len(var.shape) == 3:
            value = var[it, :, -1]  # surface value
        elif len(var.shape) == 2:
            value = var[it, :]
    elif isinstance(layer, int):  # layer number
        if len(var.shape) == 3:
            value = var[it, :, layer]
        elif len(var.shape) == 2:
            value = var[it, :]
    else:
        raise ValueError(f"Invalid layer value {layer}")

    return value


def plot_variable(
    gd, value, caxis, xlim=None, ylim=None,
    title_str=None, output_filename='output.png'
):
    """
    Plot the variable data.
    """
    plt.figure(figsize=(7, 7))
    gd.plot_grid(fmt=1, value=value, clim=caxis, levels=31, cmap='jet', xlim=xlim, ylim=ylim)
    plt.gca().set_aspect('equal', 'box')
    plt.title(title_str)
    if output_filename is not None:
        plt.savefig(output_filename, dpi=400)

    plt.show()


def calculate_disturbance(gd, value):
    """
    Calculate the disturbance variable (if needed).
    """
    disturbance = copy.deepcopy(value)
    land = gd.dp < 0.0
    disturbance[land] = value[land] + gd.dp[land]
    disturbance[~(disturbance > 0)] = -9999  # Set negative disturbance to -9999
    gd.dp = disturbance

    return disturbance


def visualize_schism_data(
    gd_fname, fnames, var_name, dry_mask_var=None, layer=-1,
    caxis=None, xlim=None, ylim=None, output_filename='output.png'
):
    """
    Main function to load, process, and visualize the Schism data.
    """
    if var_name == "disturbance":
        var_name = "elevation"
        processed_var_name = "disturbance"
    else:
        processed_var_name = var_name

    # Read data
    gd, var, dry_mask, kbp = read_schism_data(gd_fname, fnames, var_name, dry_mask_var=dry_mask_var, layer=layer)

    # Process variable
    value = process_variable(var, dry_mask, layer=layer, kbp=kbp)

    if processed_var_name == "disturbance":
        value = calculate_disturbance(gd, value)

    # Plot the variable
    plot_variable(gd, value, caxis, xlim, ylim, output_filename)


if __name__ == "__main__":
    # Configuration: Set filenames, variable names, and plotting parameters
    gd_fname = '/sciclone/schism10/feiye/STOFS3D-v8/R15c_v7/hgrid.gr3'
    var_name = "elevation"  # Variable name in the ncfile, e.g., "salt_surface" or "zeta"
    layer = 's'
    dry_mask_var = "dryFlagNode"
    fnames = ['/sciclone/schism10/feiye/STOFS3D-v8/R15d1_v7/outputs/out2d_29.nc']
    caxis = default_caxis[var_name]  # [-2, 2]
    xlim = None  # [-77, -75]
    ylim = None  # [37, 40]
    output_filename = None

    # Visualize the data
    visualize_schism_data(gd_fname, fnames, var_name, dry_mask_var, layer, caxis, xlim, ylim, output_filename)
    print('Done')
