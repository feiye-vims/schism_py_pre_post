import pandas as pd
import pickle
import os
import numpy as np
from plot_coastal_act import get_hindcast_elev
import matplotlib.pyplot as plt
from schism_py_pre_post.Download.download_usgs_with_api import download_stations, usgs_var_dict, chunks
from schism_py_pre_post.Shared_modules.generic import BinA
import pytz
from datetime import datetime


Cache_folder = os.path.realpath(os.path.expanduser('~/schism10/Cache/'))
feet2meters = 0.3048

def get_usgs_elev(station_ids=None, start_date='2021-05-01', end_date='2021-06-01'):
    downloaded_data = download_stations(
        param_id=usgs_var_dict['gauge height']['id'],
        var='guage height', station_ids=station_ids,
        datelist=pd.date_range(start=start_date, end=end_date)
    )
        
    return downloaded_data

model_start_day_str = '2021-05-01 00:00:00'
plot_start_day_str = '2021-05-01 00:00:00'
plot_end_day_str = '2021-06-01 00:00:00'
station_bp_file = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Stations/USGS_station.bp'
elev_out_files = {'RUN24a': '/sciclone/schism10/feiye/STOFS3D-v5/Outputs/O24a/fort.18',
                  'RUN01b': '/sciclone/schism10/feiye/STOFS3D-v5/Outputs/O01b_JZ/fort.18'}
mod_run_colors = ['b', 'g', 'k']

mod_run_ids = list(elev_out_files.keys())

mods = []
for id in mod_run_ids:
    elev_out_file = elev_out_files[id]
    mods.append(
        get_hindcast_elev(
            model_start_day_str=model_start_day_str,
            noaa_stations=None,
            station_in_file=station_bp_file,
            elev_out_file=elev_out_file,
            sec_per_time_unit=86400
        )
    )
time_stamps = [x.replace(tzinfo=pytz.UTC) for x in mods[0].index]
station_ids = mods[0].columns.values[:]

cache_file = f"{Cache_folder}/usgs_{plot_start_day_str.replace(' ', '_')}-{plot_end_day_str.replace(' ','_')}.pkl"
if os.path.exists(cache_file):
    with open(cache_file, 'rb') as file:
        downloaded_data = pickle.load(file)
else:
    downloaded_data = get_usgs_elev(station_ids=station_ids, start_date=plot_start_day_str, end_date=plot_end_day_str)
    with open(cache_file, 'wb') as file:
        pickle.dump(downloaded_data, file)

# arrange ids in the original order and mark missing data
downloaded_station_ids = [x.station_info['id'] for x in downloaded_data]
idx = BinA(A=station_ids, B=downloaded_station_ids)
requested_data = np.array([None] * len(station_ids))
requested_data[idx] = downloaded_data

n_subplot_row = 5
n_subplot_col = 5
plt.rcParams.update({'font.size': 10})
i_station = -1
for ichunk, chunk in enumerate(chunks(requested_data, n_subplot_col * n_subplot_row)):
    fig, ax = plt.subplots(n_subplot_col, n_subplot_row, figsize=(20, 12))
    ax = ax.ravel()
    for n, data in enumerate(chunk):
        i_station += 1
        station_id = station_ids[i_station]
        if data is not None:
            obs_date = [x.astimezone(pytz.utc) for x in data.df['date']]
            ax[n].plot(obs_date, (data.df['value'] - data.df['value'].mean()) * feet2meters, 'r.', label='obs')
            if data.station_info['id'] != station_id:
                raise Exception('Station ids for obs and mod do not match')
        for i, [mod_run_id, mod] in enumerate(zip(mod_run_ids, mods)):
            ax[n].plot(mod.index, (mod[station_id] - mod[station_id].mean()), mod_run_colors[i], label=mod_run_id)

        ax[n].title.set_text(station_id)
        ax[n].tick_params(labelrotation=20)
        ax[n].set_xlim([datetime.strptime(plot_start_day_str, "%Y-%m-%d %H:%M:%S"),
                        datetime.strptime(plot_end_day_str, "%Y-%m-%d %H:%M:%S")])
        if n == 0:
            ax[n].legend()
    for i in range(n + 1, n_subplot_col * n_subplot_row):
        ax[i].axis('off')
    plt.tight_layout(h_pad=1, w_pad=1)
    # plt.show()
    plt.savefig(f'Chunk_{ichunk}.png')
    plt.close(fig)

os.system('scp Chunk_* $cdir/srv/www/htdocs/yinglong/feiye/STOFS3D-v5/RUN24a/USGS/')
os.system('rm Chunk_*')

pass
