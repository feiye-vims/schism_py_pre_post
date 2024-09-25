"""
Functions to download water level data from NOAA CO-OPS API.
Both native api method (not dependent on any package) and
methods based on noaa_coops are provided.
"""


import os
import time
from dataclasses import dataclass
import pickle
from datetime import datetime
from json import JSONDecodeError

import numpy as np
import pandas as pd
import requests

# searvey seems to refuse to download station information if there is no data
from searvey.coops import COOPS_Station, COOPS_Query

# noaa_coops seems more tolerant
import noaa_coops  # pip install noaa-coops

from schism_py_pre_post.Utilities.util import parse_date


@dataclass
class Unify:
    '''
    A class to unify the format of downloaded elevation data
    from different packages.
    '''
    data_df_columns = {
        'searvey': {
            't': 'date_time', 'v': 'water_level', 's': 'sigma', 'f': 'flags', 'q': 'QC'
        },
        'noaa_coops': {
            't': 'date_time', 'v': 'water_level', 's': 'sigma', 'f': 'flags', 'q': 'QC'
        },
    }


def reformat_data(df, download_method='noaa_coops'):
    '''
    Reformat the downloaded data into a unified format.
    '''
    df.reset_index(inplace=True)  # move t to a column
    df = df.rename(columns=Unify.data_df_columns[download_method])
    df.set_index('date_time', inplace=True)

    return df


# -----------------------------------------------------------------------
# -------------custom functions, not dependent on other packages---------
# -----------------------------------------------------------------------
def read_json_via_api(**params):
    '''
    download json using api:

    Example url:
    https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?
    begin_date=20240301&end_date=20240301&station=8725520&
    product=water_level&datum=NAVD&time_zone=gmt&units=metric&
    application=DataAPI_Sample&format=json

    Example params:
    params = {
        "begin_date": "20240301",
        "end_date": "20240301",
        "station": "8725520",
        "product": "water_level",
        "datum": "NAVD",
        "time_zone": "gmt",
        "units": "metric",
        "application": "DataAPI_Sample",
        "format": "json"
    }
    '''
    api_url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

    # Make a GET request to the API
    # try a few times
    max_attempts = 3
    data = None
    for i in range(max_attempts):
        print(f"Attempt {i+1}")
        response = requests.get(api_url, params=params, timeout=20)

        if response is not None:
            data = response.json()

            if 'error' in data.keys():
                print(f"Error: {data['error']['message']}")
                break

            if 'data' not in data.keys() or len(data['data']) == 0:
                print("No data available for the specified parameters,"
                      f" attempts remaining: {max_attempts - i - 1}")
                time.sleep(1)
                continue
            else:
                return data  # found data, return

        else:  # retry if status code is not 200
            time.sleep(1)  # Sleep for 1 second to avoid hitting the API too frequently

    return None


def get_coops_water_level(
    begin_date='20240301', end_date='20240301',
    station='8725520', datum='NAVD', time_zone='gmt', units='metric'
):
    '''
    Get water level data from NOAA CO-OPS API at a specific station
    '''
    data = read_json_via_api(
        begin_date=begin_date, end_date=end_date, station=station,
        product="water_level", datum=datum, time_zone=time_zone, units=units,
        application="DataAPI_Sample", format="json"
    )
    if data is None:
        return None

    # Extract data points
    data_points = data['data']

    # Create lists to store datetime and water level values
    datetimes = []
    water_levels = []

    # Iterate through data points and extract datetime and water level values
    for point in data_points:
        if point['v'] == '' or point['t'] == '':
            continue  # skip if either water level or datetime is missing
        datetimes.append(point['t'])
        water_levels.append(float(point['v']))  # Convert water level to float

    # Create DataFrame
    df = pd.DataFrame({'date_time': datetimes, 'water_level': water_levels})
    df['date_time'] = pd.to_datetime(df['date_time'])
    df.set_index('date_time', inplace=True)

    return df


def sample_get_coops_water_level():
    '''
    Sample function to demonstrate how to use get_coops_water_level function
    '''
    my_df = get_coops_water_level(
        begin_date='20240312', end_date='20240326', station='8548989',
        datum='NAVD', time_zone='gmt', units='metric')
    print(my_df.head())
    my_df.to_csv('water_level_data.csv')


# -----------------------------------------------------------------------
# -------------functions dependent on other packages---------------------
# -----------------------------------------------------------------------

def get_coops_station_info(station_ids: list, method='noaa_coops'):
    '''
    Get station information from NOAA CO-OPS API.
    Wrap around noaa_coops.Station and searvey
    to handle JSONDecodeError, which occurs when
    the station does not exist or not available
    at the time of the request.

    Inputs:
    - station_ids: list of int or string

    Two methods are available:
    - 'noaa_coops': uses noaa_coops package, which takes either int or string as station id
    - 'searvey': uses searvey package, only takes int as station id
    The default method is 'noaa_coops', which is faster.

    The returned object is a Station object from either noaa_coops or searvey.
    But the following attributes are guaranteed to be available:
    - id (string)
    - name (string)
    - lon (float)
    - lat (float)

    noaa.coops.Station also has the following attributes:
    - datums
    etc.
    '''

    if method == 'noaa_coops':
        retrieve_method = noaa_coops.Station
    elif method == 'searvey':
        retrieve_method = COOPS_Station
    else:
        raise ValueError(f"Unknown method {method}. Use 'noaa_coops' or 'searvey'.")

    stations_info = []
    for station_id in station_ids:
        try:
            station_info = retrieve_method(int(station_id))
        except JSONDecodeError:
            print(f'Failed to retrieve station {station_id}')
            station_info = None

        # unify the attribute names for the two methods
        if method == 'noaa_coops':
            station_info.lon = station_info.lat_lon['lon']
            station_info.lat = station_info.lat_lon['lat']
        elif method == 'searvey':
            station_info.lon = station_info.location.coords[0][0]
            station_info.lat = station_info.location.coords[0][1]

        stations_info.append(station_info)

    return stations_info


def sample_get_coops_station_info():
    '''
    Sample function to demonstrate how to use get_coops_station_info function
    '''
    stations_info = get_coops_station_info(['8779770', 8725520], method='noaa_coops')
    for station_info in stations_info:
        if station_info is None:
            print('Station info not available')
        else:
            print(station_info.id)
            print(station_info.name)
            print(station_info.lon)
            print(station_info.lat)


def get_coops_elev(
    begin_time: datetime, end_time: datetime, noaa_stations=None,
    retrieve_method='noaa_coops', default_datum='NAVD',
    cache_folder='/sciclone/schism10/feiye/Cache/'
):

    '''
    Download elevation data for a coops station within a specific time range.
    Inputs:
    - begin_time: datetime object, but string can also be parsed
    - end_time: datetime object, but string can also be parsed
    - noaa_stations: list of station ids (string or int)
    - retrieve_method: 'noaa_coops', 'searvey', or 'native'
    - default_datum: 'NAVD' or 'MSL'
    - cache_folder: folder to save cache files, set to None to disable cache

    Outputs:
    - noaa_df_list: list of pandas DataFrames
        The data is reformatted into a pandas DataFrame to unify the format from
        different downloading methods. The columns are:
        'date_time', 'water_level', 'sigma', 'flags', 'QC'
    - datum_list: list of strings, the eventual datum used for each station
        (some requested datum may not be available, so a fallback datum is used)
    - station_data_list: list of station information objects

    "hourly_height" is used as the product for the data retrieval.
    '''

    # if not receiving datetime objects, try to parse date strings
    begin_time, _ = parse_date(begin_time)
    end_time, _ = parse_date(end_time)
    begin_time_str = begin_time.strftime('%Y-%m-%d %H:%M:%S')
    end_time_str = end_time.strftime('%Y-%m-%d %H:%M:%S')

    # lists to be returned
    noaa_df_list = []
    datum_list = []
    station_data_list = []

    for i, st in enumerate(noaa_stations):
        print(f'Processing Station {i+1} of {len(noaa_stations)}: {st}')

        # save each station into a separate cache file
        if cache_folder is not None:
            cache_filename = (
                f"{cache_folder}/COOPS_{begin_time_str.replace(' ','_')}_"
                f"{end_time_str.replace(' ','_')}_{st}_requested_{default_datum}.pkl"
            )
        else:
            cache_filename = None

        cache_success = False
        if os.path.exists(cache_filename):
            with open(cache_filename, 'rb') as f:  # Python 3: open(..., 'rb')
                this_noaa_df, this_datum, station_data = pickle.load(f)
                print(f'Existing obs data read from {cache_filename}')
                cache_success = True

        if not cache_success:
            print('Failed to read from cache, regenerating ...')

            for ntry in range(1, 6):  # try a few times to get station info
                try:
                    if retrieve_method in ['native', 'noaa_coops']:
                        station_data = noaa_coops.Station(st)
                        break  # Success, exit loop
                    elif retrieve_method == 'searvey':  # searvey, convert into the same format as noaa_coops
                        station_data = COOPS_Station(int(st))
                        lon_lat = np.squeeze(np.array(station_data.location.coords))
                        setattr(station_data, 'lat_lon', {'lat': lon_lat[-1], 'lon': lon_lat[0]})
                        break  # Success, exit loop
                    else:
                        raise ValueError(f"retrieve_method '{retrieve_method}' not supported")
                except (ConnectionError, TimeoutError) as e:
                    print(f"Network-related error for {st} on try {ntry}: {e}. Retrying...")
                except ValueError as e:
                    print(f"Value error: {e}")
                    raise  # Re-raise if it's a ValueError as it's likely not retryable

                time.sleep(2 ** ntry)  # Exponential backoff

            this_noaa_df = None
            this_datum = None
            for datum in [default_datum, 'MSL']:  # fall back to MSL, which is generally available
                print(f"Trying datum {datum}")
                try:
                    if retrieve_method == 'noaa_coops':
                        this_noaa_df = station_data.get_data(
                            begin_date=begin_time.strftime('%Y%m%d'),
                            end_date=end_time.strftime('%Y%m%d'),
                            product="hourly_height", datum=datum, units="metric", time_zone="gmt"
                        )
                        this_noaa_df = reformat_data(this_noaa_df, download_method='noaa_coops')
                    elif retrieve_method == 'searvey':
                        this_noaa_df = COOPS_Query(
                            station=st, start_date=begin_time_str, end_date=end_time_str,
                            datum=datum, product='hourly_height',
                            units='metric', time_zone='gmt'
                        ).data
                        if this_noaa_df.empty:
                            print(f"Searvey failed to get data for {st} with datum {datum}")
                            this_noaa_df = None
                        else:
                            this_noaa_df = reformat_data(this_noaa_df, download_method='searvey')
                    elif retrieve_method == 'native':  # native method from this package
                        this_noaa_df = get_coops_water_level(
                            begin_date=begin_time.strftime('%Y%m%d'), end_date=end_time.strftime('%Y%m%d'),
                            datum=datum, station=st
                        )
                    else:
                        raise ValueError(f'retrieve_method {retrieve_method} not supported')
                except (noaa_coops.station.COOPSAPIError) as e:
                    print(f"Failed to get data for {st} with datum {datum}: {e}")

                # clip data to the desired time range
                if this_noaa_df is not None:
                    this_noaa_df = this_noaa_df[begin_time:end_time]
                    if this_noaa_df.empty:
                        print(f"Some data was found for {st} with datum {datum}, but not within the desired time range")
                        this_noaa_df = None

                if this_noaa_df is not None:
                    this_datum = datum  # record actual datum retrieved
                    with open(cache_filename, 'wb') as f:
                        pickle.dump([this_noaa_df, this_datum, station_data], f)
                    break  # found data, break

        noaa_df_list.append(this_noaa_df)
        datum_list.append(this_datum)
        station_data_list.append(station_data)

    return [noaa_df_list, datum_list, station_data_list]


if __name__ == '__main__':
    sample_get_coops_station_info()
    print('Done')
