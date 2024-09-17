'''
This script is used to find the difference between two station files
'''

import pandas as pd
import numpy as np


PRODUCT1 = 'CERA'
PRODUCT2 = 'STOFS3D'

st1 = pd.read_csv(
    '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/CERA/STOFS3D_Atlantic_CERA_stations.txt',
    dtype={'stationid': str})
st2 = pd.read_csv(
    '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/CERA/station.in',
    dtype={'stationid': str})

# join the two dataframes by stationid
st1.set_index('stationid', inplace=True)
st2.set_index('stationid', inplace=True)
st = st1.join(st2, lsuffix=f'_{PRODUCT1}', rsuffix=f'_{PRODUCT2}', how='inner')

print(f'Number of stations in st1: {len(st1)}')
print(f'Number of stations in st2: {len(st2)}')
print(f'Number of stations in both: {len(st)}')

# check longitude and latitude difference
st['diff_lon'] = st[f'lon_{PRODUCT1}'] - st[f'lon_{PRODUCT2}']
st['diff_lat'] = st[f'lat_{PRODUCT1}'] - st[f'lat_{PRODUCT2}']

diff_idx = np.logical_or(abs(st['diff_lon'].values) > 1e-6, abs(st['diff_lat'].values) > 1e-6)
st_diff = st[diff_idx]

st_diff.to_csv('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/CERA/station_diff.txt', sep=',')

print('Done!')
