'''
Plot NWM time series at line segments (FIDs) of the hydrofabric
and compare with USGS time series at the corresponding USGS stations.
'''

import xarray as xr
from glob import glob

import numpy as np
import pandas as pd
from datetime import datetime
import pytz
from matplotlib import pyplot as plt

from schism_py_pre_post.Download.download_usgs_with_api import download_stations, usgs_var_dict


def b_in_a(a=None, b=None):
    '''
    Given two arrays A and B, return the indices of A's elements in B.
    '''

    if a is None and b is None:
        print('Demonstration of b_in_a')
        a = np.array([3, 5, 7, 1, 9, 8, 6, 6])
        b = np.array([3, 1, 5, 8, 6])
        print(f'A: {a}')
    else:
        a = np.array(a)
        b = np.array(b)

    index = np.argsort(a)
    sorted_a = a[index]
    sorted_index = np.searchsorted(sorted_a, b)

    bindex = np.take(index, sorted_index, mode="raise")
    mask = a[bindex] != b

    result = np.ma.array(bindex, mask=mask)
    return result


def get_nwm_var(var_str="streamflow", nwm_files=None, fids=None):
    '''
    Get a variable from NWM
    :param var_str: str, variable name
    :param nwm_dir: list, NWM output files, each file contains a time step
    :param fids: np array or list, feature IDs of interest

    :return: var: np.ndarray, shape=(n_files, n_fids), variable values
    '''
    if nwm_files is None or len(nwm_files) == 0:
        raise ValueError('nwm_files is not provided')
    if fids is None or len(fids) == 0:
        raise ValueError('fids is not provided')

    fids = np.array(fids)

    river_heads_fid = np.array(list(fids)).astype(int)
    with xr.open_dataset(nwm_files[0]) as ds:
        fid = ds['feature_id'].values
        river_heads_idx = b_in_a(fid, river_heads_fid)

    var = np.zeros((len(nwm_files), len(river_heads_fid)))
    time_stamps = np.zeros((len(nwm_files), ), dtype=object)
    for i, file in enumerate(nwm_files):
        print(file)
        ds = xr.open_dataset(file)
        var[i, :] = ds[var_str].values[river_heads_idx]
        if len(ds['time']) != 1:
            raise ValueError('Multiple time steps in a file')
        time_stamps[i] = ds['time'].values[0]

        # deal with nan
        if np.any(np.isnan(var[i, :])):
            print(f'Warning: nan found in {var_str} at requested FIDs in {file} \
                   at time {time_stamps[i]} \
                   at {len(np.where(np.isnan(var[i, :]))[0])} locations')
            var[np.isnan(var)] = 0.0
            print('Warning: nan replaced with 0.0')

    return var, time_stamps


def view_nwm():
    '''View time-series of NWM variables at line segments (FIDs) of the hydrofabric'''
    # usgs2fid_dict = {
    #     '02492511': ['15707937', '15708315'],
    #     '02492600': ['15708755'],
    #     '07375050': ['18931936'],
    #     '07375170': ['18928090'],
    #     '07375500': ['18975531'],
    #     '07376000': ['20090368'],
    #     '07378050': ['18988626'],
    #     '07378500': ['18990204'],
    # }
    # inputs
    fids = ['15720717', '15707937', '15721007', '15707927', '15707961']
    fids = ['15148368', '21898499', '21898517']
    fids = ['15148144']
    fids = ['19406836']  # , '19406834', '19406822', '19406814']  # Atchafalaya

    fid2USGS_dict = {
        '19406836': '07381490',
        # '19406822': '07381490',
    }
    fids = list(fid2USGS_dict.keys())

    fids = ['18991160', '18991158', '18991148', '18991146', '18991150']
    fid2USGS_dict = {
        '18991160': '07379075',
        '18991158': '07379075',
        '18991148': '07379075',
        '18991146': '07379075',
        '18991150': '07379075',
    }
    # fids = ['18991608']

    fids = ['19269176']  # Mississippi River
    fid2USGS_dict = {
        '19269176': '07374000',
    }
    # -------------------- plot NWM time series --------------------
    # get nwm_files
    # nwm_files = sorted(glob('/sciclone/schism10/feiye/STOFS3D-v8/I09/Source_sink/'
    #                         'original_source_sink/20240305/nwm.t00z.medium_range.channel_rt_1.*.nc'))
    nwm_files = sorted(glob(
        '/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/CHRTOUT/for_2019_hindcast/201906090100.CHRTOUT_DOMAIN1'
        # '/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/CHRTOUT/2019/2019*CHRTOUT*'
    ))
    # end inputs

    nwm, time_stamps = get_nwm_var(var_str="streamflow", nwm_files=nwm_files, fids=fids)

    # plot NWM time series
    _, ax = plt.subplots()
    for i, fid in enumerate(fids):
        ax.plot(time_stamps, nwm[:, i], label=f'NWM fid: {fid}')
        ax.legend()

    # -------------------- plot USGS at the corresponding FIDs --------------------

    # pad one day before and after the start and end time
    padded_start_time = time_stamps[0] - pd.Timedelta('1 day')
    padded_end_time = time_stamps[-1] + pd.Timedelta('1 day')
    usgs_var = 'gauge height'
    additional_scale = 10  # to be plotted in the same scale as NWM flow
    usgs_data = download_stations(
        param_id=usgs_var_dict[usgs_var]['id'],
        station_ids=list(fid2USGS_dict.values()),
        datelist=pd.date_range(start=padded_start_time, end=padded_end_time),
    )
    usgs_data_dict = {}
    for i, data in enumerate(usgs_data):
        usgs_data_dict[data.station_info['id']] = data.df

    for i, usgs_id in enumerate(np.unique(list(fid2USGS_dict.values()))):
        usgs_df = usgs_data_dict[usgs_id]
        usgs_time = pd.to_datetime(usgs_df["date"], utc=True)
        ax.plot(usgs_time, usgs_df['value'] * usgs_var_dict[usgs_var]['unit_conv'] * additional_scale, '--r', label=f'USGS: {usgs_id}')
        ax.legend()
    
    plt.rcParams.update({'font.size': 16}) 
    plt.grid()
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.xlim(datetime(2024, 3, 12), datetime(2024, 4, 6))
    plt.show()
    print('Done')


if __name__ == '__main__':
    view_nwm()
    print('Done')
