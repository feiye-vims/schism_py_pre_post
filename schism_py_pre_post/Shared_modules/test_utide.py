"""Testing the utide package with SCHISM time series data."""

# import numpy as np
from pylib_experimental.schism_file import TimeHistory
from utide import solve


my_th = TimeHistory.from_file(
    '/sciclone/schism10/feiye/STOFS3D-v8/O04c/elevation.USGS_station_LA_repositioned.dat',
    th_unit='days', start_time_str='2024-03-05 00:00:00')

ha_results = [None] * my_th.n_station
for i in range(25):
    ha_results[i] = solve(t=my_th.datetime, u=my_th.data[:, i], lat=30)

for i in range(25):
    # idx = np.argmax(ha_results[i].A)
    # print(f'{my_th.stations[i]}: {ha_results[i].name[idx]} amplitude = {ha_results[i].A[idx]} m')

    idx = ha_results[i].name == 'O1'
    print(f'{my_th.stations[i]}: O1 amplitude = {ha_results[i].A[idx]} m, '
          f'phase = {ha_results[i].g[idx]} deg')

print('Done')
