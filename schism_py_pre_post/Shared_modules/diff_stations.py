'''
This script is used to find the difference between two station files
'''

from pathlib import Path
import pickle
from json import JSONDecodeError
import pandas as pd
import numpy as np
from searvey.coops import COOPS_Query, COOPS_Station
import noaa_coops


def diff_cera():
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


def get_original_stations(station_file: Path = None):
    if station_file is None:
        station_file = Path('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/CERA/station.in')

    st = pd.read_csv(station_file, dtype={'stationid': str})
    stnames = list(st['stationid'].values.astype(int))[:164]
    
    st_list = []
    for stname in stnames:
        print(f'retrieving station {stname}')
        try:
            stinfo = noaa_coops.Station(stname)
        except JSONDecodeError:
            print(f'Failed to retrieve station {stname}')
            stinfo = None

        st_list.append(stinfo)

    with open('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/CERA/station_info.pkl', 'wb') as file:
        pickle.dump(st_list, file)


if __name__ == '__main__':
    get_original_stations()