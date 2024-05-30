#!/usr/bin/env python
"""Replace a source"""
import numpy
import pandas as pd
import pytz
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

# pip install climata
from schism_py_pre_post.Download.download_usgs_with_api import download_single_station
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Prop import Prop


vsource_file = '/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Vsource/NOAA/vsource.th'
vsource_file = '/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Vsource/VIMS/vsource.th'
start_time_str = '2022-08-15 00:00:00'  # add one day before

time_fmt = '%Y-%m-%d %H:%M:%S'
source_eles = Prop('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v4/Scale_vsource/source_sink.in', 2).ip_group[0]  # first group is source
my_sources = TimeHistory(vsource_file, start_time_str=start_time_str)
end_time_str = my_sources.df.iloc[-1, 0].strftime(time_fmt)
prev_scales = [1.0]

Missi_idx = my_sources.get_max_station()
Missi_ele_id = source_eles[Missi_idx]

for ele, usgs_station_id, prev_scale in zip([Missi_ele_id], ['07374000'], prev_scales):
    # ----------------------------------------------------------------------------------
    # -----------------replace a source-------------------------------------------------------
    # ----------------------------------------------------------------------------------
    # obs
    # pad start and end time
    obs_start_time_str = (datetime.strptime(start_time_str, time_fmt) - timedelta(days=1)).strftime(time_fmt)
    obs_end_time_str = (datetime.strptime(end_time_str, time_fmt) + timedelta(days=1)).strftime(time_fmt)
    data = download_single_station(station_id=usgs_station_id, param_id='00060',
                                   var='streamflow',
                                   datelist=pd.date_range(start=obs_start_time_str, end=obs_end_time_str))
    obs_df = data[0].df
    # convert time to utc
    obs_df['date'] = [x.tz_convert(pytz.utc) for x in obs_df['date']]
    obs_time = (obs_df['date'] - pd.Timestamp(datetime.strptime(start_time_str, time_fmt), tz='UTC'))
    obs_sec = [x.total_seconds() for x in obs_time]
    unit_conv = 0.028316846592
    obs_val = obs_df['value'].values * unit_conv

    # ----------change the vsources at specified eles-----------
    eles_replaced = [ele]  # eles_added = ['2096361']
    for ele_replaced in eles_replaced:
        idx = numpy.where(source_eles == int(ele_replaced))
        if len(idx[0]) != 1 and idx[0][0] < 0:
            raise Exception('idx error')
        else:
            idx = idx[0][0]
        # my_sources.plot_ts(idx, 'source', 1)

        # from a time series
        vs_obs = numpy.interp(
            numpy.array(my_sources.time, dtype=float),
            numpy.array(obs_sec, dtype=float),
            obs_val.astype(numpy.float64)
        )

        vs_obs[vs_obs < 0.01] = 0.01
        scale = vs_obs.mean() / (my_sources.df.iloc[:, idx + 1].mean() / prev_scale)

        plt.plot(my_sources.datetime, my_sources.df.iloc[:, idx + 1], label='NWM')
        my_sources.df.iloc[:, idx + 1] = vs_obs * 1.0
        plt.plot(my_sources.datetime, my_sources.df.iloc[:, idx + 1], label='Obs')
        plt.legend()
        plt.show()

        print(f'previous scale: {prev_scale}')
        print(f'source_ele_idx: {idx}; scale: {scale}')

        # from a simpler formulae
        # vs_added = my_sources.data[:, idx]
        # vs_added[vs_added < 0.01] = 0.01
        # my_sources.data[:, idx] = vs_added

    my_sources.writer(f'{my_sources.source_file}.modified')
