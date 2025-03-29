import os
import pickle
from datetime import datetime

import pytz
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from schism_py_pre_post.Plot.plot_elev import get_hindcast_elev
from schism_py_pre_post.Download.download_usgs_with_api import \
    download_stations, usgs_var_dict, chunks
from schism_py_pre_post.Shared_modules.generic import BinA


FEET2METERS = 0.3048


def get_usgs_data(
    station_ids=None, var='gauge height',
    start_date='2021-05-01', end_date='2021-06-01',
    cache_ident_name=None
):
    '''
    Get USGS data for a list of stations, with cache handling
    '''
    # handle cache: cache file is saved per station group specified in a bp file.
    if cache_ident_name is not None:
        cache_folder = os.path.realpath(os.path.expanduser('~/schism10/Cache/'))
        cache_file = (f"{cache_folder}/usgs_{var}_{cache_ident_name}_"
                      f"{start_date.replace(' ', '_')}-{end_date.replace(' ','_')}.pkl")

        if os.path.exists(cache_file):
            with open(cache_file, 'rb') as file:
                downloaded_data = pickle.load(file)

    if cache_ident_name is None or not os.path.exists(cache_file):
        downloaded_data = download_stations(
            param_id=usgs_var_dict[var]['id'],
            station_ids=station_ids,
            datelist=pd.date_range(start=start_date, end=end_date)
        )

    if cache_ident_name is not None:
        with open(cache_file, 'wb') as file:
            pickle.dump(downloaded_data, file)

    # Since the downloader may change the station order in station.bp,
    # arrange ids in the original order and mark missing data
    downloaded_station_ids = [x.station_info['id'] for x in downloaded_data]
    idx = BinA(A=station_ids, B=downloaded_station_ids)
    requested_data = np.array([None] * len(station_ids))
    requested_data[idx] = downloaded_data

    return requested_data


def get_model_elev(elev_out_files, model_start_day_str, station_bp_file, sec_per_time_unit):
    """
    Get model results for multiple runs
    """
    mod_run_ids = list(elev_out_files.keys())

    mods = []
    for run_id in mod_run_ids:
        elev_out_file = elev_out_files[run_id]
        mods.append(
            get_hindcast_elev(
                model_start_day_str=model_start_day_str,
                noaa_stations=None,
                station_in_file=station_bp_file,
                elev_out_file=elev_out_file,
                sec_per_time_unit=sec_per_time_unit,
            )
        )

    return mod_run_ids, mods


def plot_elev_no_stats(
    mods, mod_run_ids, requested_obs_data,
    plot_start_day_str, plot_end_day_str,
    output_dir='./', mod_run_colors=None
):
    """
    Plot model results and observations without statistics
    """
    if mod_run_colors is None:
        mod_run_colors = ['c', 'k', 'm', 'y', 'g', 'b']

    station_ids = mods[0].columns.values[:]  # a list of station ids
    # time_stamps = [x.replace(tzinfo=pytz.UTC) for x in mods[0].index]

    # plot
    n_subplot_row = 7
    n_subplot_col = 4
    n_chunk = 25
    plt.rcParams.update({'font.size': 20})
    i_station = -1
    for ichunk, chunk in enumerate(chunks(requested_obs_data, n_chunk)):
        fig, ax = plt.subplots(n_subplot_row, n_subplot_col, figsize=(60, 30))
        ax = ax.ravel()
        for n, data in enumerate(chunk):
            i_station += 1
            station_id = station_ids[i_station]
            if data is not None:
                obs_date = [x.astimezone(pytz.utc) for x in data.df['date']]
                obs_y = (data.df['value'] - data.df['value'].mean()) * FEET2METERS
                ax[n].plot(obs_date, obs_y, 'r.', label='obs')
                if data.station_info['id'] != station_id:
                    raise ValueError('Station ids for obs and mod do not match')
            for i, [mod_run_id, mod] in enumerate(zip(mod_run_ids, mods)):
                mod_y = mod[station_id] - mod[station_id].mean()
                ax[n].plot(mod.index, mod_y, mod_run_colors[i], label=mod_run_id)

            ax[n].title.set_text(station_id)
            ax[n].tick_params(labelrotation=20)
            ax[n].set_xlim([datetime.strptime(plot_start_day_str, "%Y-%m-%d %H:%M:%S"),
                            datetime.strptime(plot_end_day_str, "%Y-%m-%d %H:%M:%S")])
            if data is not None and not data.df.empty:  # set ylim to obs's range
                ax[n].set_ylim(obs_y.min() - 0.5, obs_y.max() + 0.5)
            if n == 0:
                ax[n].legend()
            ax[n].grid()
            # for i in range(n + 1, n_subplot_col * n_subplot_row):
            #     ax[i].axis('off')
        plt.tight_layout(h_pad=1, w_pad=1)
        # plt.show()
        plt.savefig(f'Chunk_{ichunk}.png')
        plt.close(fig)

    os.makedirs(output_dir, exist_ok=True)
    os.system(f'scp Chunk_* {output_dir}')
    os.system('rm Chunk_*')


def plot_usgs(
    station_bp_file=None, model_start_day_str='2021-05-01 00:00:00',
    plot_start_day_str='2021-05-01 00:00:00', plot_end_day_str='2021-06-01 00:00:00',
    output_dir=None, elev_out_files: dict = None, sec_per_time_unit=86400,
):
    '''
    Plot USGS data and model results

    mulitple model runs are supported, e.g.,
    elev_out_files = {
        'RUN24a': '/sciclone/schism10/feiye/STOFS3D-v5/Outputs/O24a/fort.18',
        'RUN01b': '/sciclone/schism10/feiye/STOFS3D-v5/Outputs/O01b_JZ/fort.18'
    },
    '''

    print('Gathering model results')
    mod_run_ids, mods = get_model_elev(
        elev_out_files=elev_out_files,
        model_start_day_str=model_start_day_str,
        station_bp_file=station_bp_file,
        sec_per_time_unit=sec_per_time_unit,
    )

    print('Gathering observation')
    obs_data = get_usgs_data(
        station_ids=mods[0].columns.values[:],  # a list of station ids
        start_date=plot_start_day_str,
        end_date=plot_end_day_str,
        cache_ident_name=os.path.basename(station_bp_file)
    )

    print('Plotting')
    plot_elev_no_stats(
        mods=mods,
        mod_run_ids=mod_run_ids,
        requested_obs_data=obs_data,
        plot_start_day_str=plot_start_day_str,
        plot_end_day_str=plot_end_day_str,
        output_dir=output_dir,
    )


scenarios_dict = {
    'v8_March_reforecast_LA': {
        'station_bp_file': '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA.bp',
        'model_start_day_str': '2024-03-05 00:00:00',
        'plot_start_day_str': '2024-03-10 00:00:00',
        'plot_end_day_str': '2024-04-10 00:00:00',
    },
    'v8_March_reforecast_LA_nontidal': {
        'station_bp_file': '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal_v43.bp',
        'model_start_day_str': '2024-03-05 00:00:00',
        'plot_start_day_str': '2024-03-10 00:00:00',
        'plot_end_day_str': '2024-04-10 00:00:00',
    },
    'Ida_LA': {
        'station_bp_file': '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA.bp',
        'model_start_day_str': '2021-08-01 00:00:00',
        'plot_start_day_str': '2021-08-01 00:00:00',
        'plot_end_day_str': '2021-09-10 00:00:00',
    },
    'v8': {
        'station_bp_file': '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_repositioned_v49.bp',
        'model_start_day_str': '2024-03-05 00:00:00',
        'plot_start_day_str': '2024-03-10 00:00:00',
        'plot_end_day_str': '2024-04-10 00:00:00',
    },
    '2018_hindcast': {
        'station_bp_file': '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_repositioned_v49.bp',
        'model_start_day_str': '2017-12-01 00:00:00',
        'plot_start_day_str': '2018-07-01 00:00:00',
        'plot_end_day_str': '2018-08-12 00:00:00',
    },
}


def viz_usgs():
    '''quick visualization of usgs data'''

    usgs_stations = ['01021060', '01022840', '011058837', '02299734']
    usgs = get_usgs_data(
        station_ids=usgs_stations, start_date='2024-03-05', end_date='2024-04-10',
        cache_ident_name=None, var='gauge height')

    # plot time series in nx1 subplot
    cubic_feet2cubic_meters = 0.0283168
    _, ax = plt.subplots(len(usgs_stations)+1, 1, figsize=(10, 5))
    for i, station in enumerate(usgs_stations):
        ax[i].plot(
            usgs[i].df['date'], usgs[i].df['value']*cubic_feet2cubic_meters,
            'r.', label='USGS'
        )
        ax[i].set_title(station)
        ax[i].legend()
    plt.show()


if __name__ == "__main__":

    scenario = scenarios_dict['2018_hindcast']  # see scenarios_dict for options

    plot_usgs(
        station_bp_file=scenario['station_bp_file'],
        model_start_day_str=scenario['model_start_day_str'],
        sec_per_time_unit=86400,
        # elev_out_files={
        #     'R15a_v7': '/sciclone/schism10/feiye/STOFS3D-v8/O15a_v7/elevation.USGS_station_LA_repositioned_nontidal.dat',
        #     'R09c': '/sciclone/schism10/feiye/STOFS3D-v8/O09c/elevation.USGS_station_LA_repositioned_nontidal.dat',
        # },
        elev_out_files={
            'RUN16_v6': '/sciclone/schism10/feiye/STOFS3D-v8/O16_v6/elevation.USGS_repositioned_v49.dat',
            'R20b': '/sciclone/schism10/feiye/STOFS3D-v8/O20b/elevation.USGS_repositioned_v50.dat',
        },
        plot_start_day_str=scenario['plot_start_day_str'],
        plot_end_day_str=scenario['plot_end_day_str'],
        output_dir='/sciclone/schism10/feiye/STOFS3D-v8/O20b/',
    )
