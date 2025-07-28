"""
Daily automated plotting script for the operation of the STOfS3D system.
"""

import os
from pathlib import Path
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone, timedelta
import boto3
from botocore import UNSIGNED
from botocore.config import Config
from pylib import schism_grid as read_schism_hgrid


oper_dir = Path('/sciclone/schism10/feiye/STOFS3D_Crontab/')


def assemble_forecast_dict(forecast_date, output_dir, forecast_version):
    """
    make a dictionary of variables to plot with their parameters

    inputs:
    - forecast_date_str: string, date in 'YYYYMMDD' format
    - output_dir: Path, directory to save the plots
    - forecast_version: string, version of the forecast
    """

    forecast_date_str = forecast_date.strftime("%Y%m%d")
    forecast_dict = {'forecast_date_str': forecast_date_str}

    if forecast_version == 'v3.1.0':
        forecast_dict['parent_key'] = f'STOFS-3D-Atl/para/VIMS_v72_stofs.v3.1.0/stofs_3d_atl.{forecast_date_str}'
        forecast_dict['hgrid_file'] = '/sciclone/schism10/feiye/STOFS3D-v8/v7.2_static_inputs_2025_04_12/hgrid.gr3'
    elif forecast_version == 'v2.1':
        forecast_dict['parent_key'] = f'STOFS-3D-Atl/stofs_3d_atl.{forecast_dict["forecast_date_str"]}'
        forecast_dict['hgrid_file'] = '/sciclone/home/feiye/s1/STOFS3D-v6/Inputs/v6.1.static_inputs_snapshot20231105/hgrid.gr3'

    output_dir = output_dir / forecast_version / forecast_date_str
    os.makedirs(output_dir, exist_ok=True)
    forecast_dict['output_dir'] = output_dir

    plot_var_dict = {
        'temp_surface': {
            'caxis': [-2, 35],
            'title': f'STOFS3D Surface Temperature on {forecast_date_str}',
            'output_filename': f'{output_dir}/stofs3d_temp_surface_{forecast_date_str}.png'
        },
        'temp_bottom': {
            'caxis': [-2, 35],
            'title': f'STOFS3D Bottom Temperature on {forecast_date_str}',
            'output_filename': f'{output_dir}/stofs3d_temp_bottom_{forecast_date_str}.png'
        },
        'salt_surface': {
            'caxis': [0, 40],
            'title': f'STOFS3D Surface Salinity on {forecast_date_str}',
            'output_filename': f'{output_dir}/stofs3d_salt_surface_{forecast_date_str}.png'
        },
        'salt_bottom': {
            'caxis': [0, 40],
            'title': f'STOFS3D Bottom Salinity on {forecast_date_str}',
            'output_filename': f'{output_dir}/stofs3d_salt_bottom_{forecast_date_str}.png'
        }
    }
    if forecast_version == 'v3.1.0':
        plot_var_dict['zeta'] = {
            'caxis': [-2, 2],
            'title': f'STOFS3D Surface Elevation on {forecast_date_str}',
            'output_filename': f'{output_dir}/stofs3d_zeta_{forecast_date_str}.png'
        }
    elif forecast_version == 'v2.1':
        plot_var_dict['elev'] = {
            'caxis': [-2, 2],
            'title': f'STOFS3D Surface Elevation on {forecast_date_str}',
            'output_filename': f'{output_dir}/stofs3d_elev_{forecast_date_str}.png'
        }
    forecast_dict['plot_var_dict'] = plot_var_dict

    return forecast_dict


def download_stofs3d_from_s3(
    forecast_date: datetime = None, parent_key: str = None,
    forecast_version: str = 'v3.1.0'
):
    """
    Download the latest STOFS3D data from S3 bucket.

    The file is expected to be in the format:
    'STOFS-3D-Atl/para/VIMS_v72_stofs.v3.1.0/stofs_3d_atl.YYYYMMDD/stofs_3d_atl.t12z.field2d_f001_012.nc'
    or
    'STOFS-3D-Atl/stofs_3d_atl.YYYYMMDD/stofs_3d_atl.t12z.field2d_f001_012.nc'

    Args:
        forecast_date (datetime): The date for which to download the data.
        parent_key (str): The S3 key prefix where the files are located.
        forecast_version (str): The version of the forecast data.
    """

    day_str = forecast_date.strftime("%Y%m%d")
    # Construct URL with today's date
    s3 = boto3.client('s3', region_name='us-east-1', config=Config(signature_version=UNSIGNED))
    bucket = 'noaa-nos-stofs3d-pds'

    # Local file path
    local_folder = f"{oper_dir}/{forecast_version}/{day_str}"
    os.makedirs(local_folder, exist_ok=True)  # Ensure the folder exists
    
    local_files = []
    for file in ['stofs_3d_atl.t12z.field2d_f001_012.nc']:
        remote_key = f'{parent_key}/{file}'
        local_file = f"{local_folder}/{file}"

        if os.path.exists(local_file):
            print(f"File {local_file} already exists. Skipping download.")
        else:
            print(f"Downloading {remote_key} from bucket {bucket} to {local_file} ...")
            s3.download_file(bucket, remote_key, local_file)
        local_files.append(local_file)

    return local_files


def get_2d_field(filename, var_name='temp_surface', it=0):
    """
    Plot a 2D field from the STOFS3D data.
    """
    ds = xr.open_dataset(filename)
    
    if var_name not in ds:
        raise ValueError(f"Variable '{var_name}' not found in the dataset.")
    
    data = ds[var_name].isel(time=it)

    return np.array(data)


def plot_2d_field(
    gd, value, caxis, xlim=None, ylim=None,
    title_str=None, output_filename='output.png'
):
    """
    Plot 2d field data on a Schism grid.
    """
    plt.figure(figsize=(7, 7))
    gd.plot_grid(fmt=1, value=value, clim=caxis, levels=31, cmap='jet', xlim=xlim, ylim=ylim)
    plt.gca().set_aspect('equal', 'box')
    plt.title(title_str)
    if output_filename is not None:
        plt.savefig(output_filename, dpi=400)

    # plt.show()


def print_date():
    """
    Print the current date in UTC.
    """
    today = datetime.now(timezone.utc)
    print("\n------------------------------")
    print(f"Running {__file__} on {today.strftime('%Y-%m-%d %H:%M:%S')} UTC")


if __name__ == "__main__":  
    print_date()

    today = datetime.now(timezone.utc)
    yesterday = (datetime.now(timezone.utc) - timedelta(days=1))

    forecast_date = yesterday
    forecast_date_str = forecast_date.strftime("%Y%m%d")

    for forecast_version in ['v3.1.0', 'v2.1']:
        print(f'Processing STOFS3D data for version {forecast_version}...')
    
        forecast_dict = assemble_forecast_dict(forecast_date, oper_dir, forecast_version)

        print(f'Dowloading STOFS3D data for {forecast_date_str}...')
        downloaded_files = download_stofs3d_from_s3(
            forecast_date, forecast_dict['parent_key'], forecast_version)

        print('Loading STOFS3D grid...')
        hgrid = read_schism_hgrid(forecast_dict['hgrid_file'])

        for var_name, params in forecast_dict['plot_var_dict'].items():
            print(f'Plotting {var_name}...')
            field2d_file = downloaded_files[0]
            data = get_2d_field(field2d_file, var_name=var_name, it=0)
            plot_2d_field(
                hgrid, data, caxis=params['caxis'],
                title_str=params['title'],
                output_filename=params['output_filename']
            )

        print(f'Done processing STOFS3D data for version {forecast_version}.\n')


