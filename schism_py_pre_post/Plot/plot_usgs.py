import pandas as pd
import pickle
import os
from plot_coastal_act import get_hindcast_elev
import matplotlib.pyplot as plt
from schism_py_pre_post.Download.download_usgs_with_api import download_stations, usgs_var_dict, chunks
from schism_py_pre_post.Shared_modules.generic import BinA
import pytz
from datetime import datetime


feet2meters = 0.3048

def get_usgs_elev(station_ids=None, start_date='2021-05-01', end_date='2021-06-01'):
    total_data = download_stations(
        param_id=usgs_var_dict['gauge height']['id'],
        var='guage height', station_ids=station_ids,
        datelist=pd.date_range(start=start_date, end=end_date)
    )
        
    return total_data

model_start_day_str = '2021-05-01 00:00:00'
plot_start_day_str = '2021-07-24 00:00:00'
plot_end_day_str = '2021-08-01 00:00:00'
station_bp_file = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Stations/USGS_station.bp'
elev_out_file = '/sciclone/schism10/feiye/STOFS3D-v5/Outputs/O24a/fort.18'

mod = get_hindcast_elev(
    model_start_day_str=model_start_day_str,
    noaa_stations=None,
    station_in_file=station_bp_file,
    elev_out_file=elev_out_file,
    sec_per_time_unit=86400
)
time_stamps = [x.replace(tzinfo=pytz.UTC) for x in mod.index]
station_ids = mod.columns.values[:]

cache_file = f"usgs_{plot_start_day_str.replace(' ', '_')}-{plot_end_day_str.replace(' ','_')}.pkl"
if os.path.exists(cache_file):
    with open(cache_file, 'rb') as file:
        total_data = pickle.load(file)
else:
    total_data = get_usgs_elev(station_ids=station_ids, start_date=plot_start_day_str, end_date=plot_end_day_str)
    with open(cache_file, 'wb') as file:
        pickle.dump(total_data, file)
# arrange ids in the original order and mark missing data
# downloaded = [x.station_info['id'] for x in total_data]
# idx = BinA(A=station_ids, B=downloaded)

n_subplot_row = 5
n_subplot_col = 5
plt.rcParams.update({'font.size': 10})
for ichunk, total_data_chunk in enumerate(chunks(total_data, n_subplot_col * n_subplot_row)):
    fig, ax = plt.subplots(n_subplot_col, n_subplot_row, figsize=(20, 12))
    ax = ax.ravel()
    for n, data in enumerate(total_data_chunk):
        obs_date = [x.astimezone(pytz.utc) for x in data.df['date']]
        ax[n].plot(obs_date, data.df['value'] * feet2meters, 'r.', label='obs')
        ax[n].plot(mod.index, mod[data.station_info['id']], 'b', label='mod')
        ax[n].title.set_text(data.station_info['id'])
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

os.system('scp Chunk_* $cdir/srv/www/htdocs/yinglong/feiye/STOFS3D-v5/RUN24a/USGS/')
os.system('rm Chunk_*')

pass
