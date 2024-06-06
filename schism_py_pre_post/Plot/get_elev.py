'''
This script is to get elevation time series from different sources,
process them into the same format utilizing pandas dataframe,
which can be used for comparison and plotting.
'''

import os
import glob
import pickle

import math
import json
import numpy as np
import pandas as pd
from scipy.interpolate import griddata  # , interp2d
from datetime import timedelta
from datetime import datetime

import noaa_coops
from searvey.coops import COOPS_Query, COOPS_Station

from schism_py_pre_post.Grid.Bpfile import Bpfile
from schism_py_pre_post.Shared_modules.test_modules import HYCOM
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Shared_modules.obs_mod_comp import obs_mod_comp
from schism_py_pre_post.Download.download_usgs_with_api import download_stations, usgs_var_dict, chunks
from schism_py_pre_post.Utilities.util import b_in_a, parse_date


def get_hycom_elev(point_xy:np.ndarray, hycom_file=None):
    '''
    Get elevation at point_xy (nx2 numpy array)
    from a HYCOM file
    '''
    st_lon = point_xy[:, 0]
    st_lat = point_xy[:, 1]

    my_hycom = HYCOM(hycom_file)
    elev = my_hycom.get_var('surf_el')
    elev[elev < -29999] = np.nan
    lon = my_hycom.lon
    lat = my_hycom.lat

    # f = interp2d(lon, lat, elev[0, :, :], kind='linear')
    # z_interp = f(st_lon, st_lat)

    X, Y = np.meshgrid(lon, lat)
    Z = elev.reshape(elev.shape[0], -1).T
    x = X.ravel()
    y = Y.ravel()
    z = Z
    z_interp = griddata(np.c_[x, y], z, (st_lon, st_lat), method='linear')

    # dealing with nan
    nan_idx = np.where(np.isnan(z_interp[:, 0]))[0]
    idx = (Z[:, 0] > -29999)
    x = X.ravel()[idx]
    y = Y.ravel()[idx]
    z = Z[idx, :]
    z_interp_nearest = griddata(np.c_[x, y], z, (st_lon, st_lat), method='nearest')
    z_interp[nan_idx] = z_interp_nearest[nan_idx]

    # fig, ax = plt.subplots(nrows=1, ncols=1)
    # ax.pcolormesh(X, Y, elev[300, :, :], vmin=-0.2, vmax=0.2, cmap='jet')
    # ax.scatter(st_lon, st_lat, 10, z_interp[:, 300], cmap='jet', vmin=-0.2, vmax=0.2)
    # plt.show()

    model_df = pd.DataFrame(data=z_interp.T, index=my_hycom.t_datetime, columns=noaa_stations)
    return model_df


def get_hindcast_elev(model_start_day_str, noaa_stations=None, station_in_file=None, elev_out_file=None, sec_per_time_unit=1, station_in_subset=None):
    '''
    Get model time series at noaa_stations,
    either from specified ids or from station_in_file
    '''
    if noaa_stations is not None:
        pass
    else:
        st_id = Bpfile(station_in_file, cols=5).make_dataframe().columns

        my_th = TimeHistory(elev_out_file, model_start_day_str, -9999, sec_per_time_unit=sec_per_time_unit)

        if station_in_subset is not None:
            my_th = my_th.export_subset(station_idx=station_in_subset)
        
        model_df = my_th.df.set_index('datetime')

        model_df.columns = st_id
    return model_df

def get_forecast_elev(plot_start_day_str, forecast_end_day_str, fcst_folder=None, station_in_file=None, i_nowcast=False):
    '''
    forecast_end_day_str corresponds to the folder name and there may be 2 days after the date in the folder name
    i_nowcast=True: take the 1st day of each 3-day forecasts;
    i_nowcast=False: take the 2nd day of each 3-day forecasts;
    '''
    st_id = Bpfile(station_in_file, cols=5).make_dataframe().columns
    if fcst_folder is None:
        fcst_folder = '/sciclone/schism10/hyu05/NOAA_NWM/oper_3D/fcst/'

    model_df = pd.DataFrame()
    model_date_range = pd.date_range(start=plot_start_day_str, end=forecast_end_day_str, freq='D')
    for i, this_day in enumerate(model_date_range):
        date_str = datetime.strftime(this_day, "%Y%m%d")
        date_str2 = datetime.strftime(this_day+timedelta(days=-1), "%Y-%m-%d %H:%M:%S")
        th_file = glob.glob(f'{fcst_folder}/{date_str}*/staout_1')
        if len(th_file) != 1:
            raise Exception(f'th_file not found in {fcst_folder}/{date_str}*: {th_file}')
        my_th = TimeHistory(th_file[0], date_str2, -9999)
        day1 = datetime.strftime(this_day - int(i_nowcast) * timedelta(days=1), "%Y-%m-%d %H:%M:%S")
        if i == len(model_date_range) - 1:  # get one more day at the end
            day2 = datetime.strftime(this_day + (int(not i_nowcast) + 1) * timedelta(days=1) - timedelta(days=0.001), "%Y-%m-%d %H:%M:%S")
        else:
            day2 = datetime.strftime(this_day + int(not i_nowcast) * timedelta(days=1) - timedelta(days=0.001), "%Y-%m-%d %H:%M:%S")
        model_df = model_df.append(my_th.df.set_index('datetime')[day1:day2])
    # ensure uniqueness
    model_df = model_df[~model_df.index.duplicated(keep='last')]

    model_df.columns = st_id

    return model_df

def get_usgs_elev(station_ids=None, start_date='2021-05-01', end_date='2021-06-01', cache_ident_name='unspecified'):
    # handle cache: cache file is saved per stations specified in a bp file.
    # The download function can do batch download for all stations, so it is not straightforward to save a per-station cache
    Cache_folder = os.path.realpath(os.path.expanduser('~/schism10/Cache/'))
    cache_file = f"{Cache_folder}/usgs_{cache_ident_name}_{start_date.replace(' ', '_')}-{end_date.replace(' ','_')}.pkl"

    if os.path.exists(cache_file):
        with open(cache_file, 'rb') as file:
            downloaded_data = pickle.load(file)
    else:
        downloaded_data = download_stations(
            param_id=usgs_var_dict['gauge height']['id'],
            station_ids=station_ids,
            datelist=pd.date_range(start=start_date, end=end_date)
        )

        with open(cache_file, 'wb') as file:
            pickle.dump(downloaded_data, file)
        
    # Since the downloader may change the station order in station.bp,
    # arrange ids in the original order and mark missing data
    downloaded_station_ids = [x.station_info['id'] for x in downloaded_data]
    idx = b_in_a(A=station_ids, B=downloaded_station_ids)
    requested_data = np.array([None] * len(station_ids))
    requested_data[idx] = downloaded_data

    obs_df = []; st_info = []
    for data in requested_data:
        data.df.columns = ['datetime', 'water_level']
        data.df.set_index('datetime', inplace=True)
        obs_df.append(data.df)
        # reformat station_info
        station_info = {'id': data.station_info['id'], 'site_name': data.station_info['name'], 'latitude': data.station_info['lat'], 'longitude': data.station_info['lon']}
        st_info.append(station_info)

    return obs_df, st_info

def get_coops_elev(
    retrieve_method='noaa_coops',  # 'searvey' or 'noaa_coops'
    plot_start_day_str=None, plot_end_day_str=None, noaa_stations=None,
    default_datum='NAVD', cache_folder='/sciclone/schism10/feiye/Cache/'
):
    '''
    Get COOPS elevation time series at noaa_stations.
    Try to read from cache first, if not found, retrieve from COOPS server.
    The output is a list of pandas dataframe, each corresponding to a station,
    because the time range may be different for different stations due to 
    data availability.
    '''

    noaa_df_list = []
    datum_list = []
    station_data_list = []
    st_info = []
    for i, st in enumerate(noaa_stations):
        print(f'Processing Station {i+1} of {len(noaa_stations)}: {st}')

        cache_filename = f"{cache_folder}/COOPS_{plot_start_day_str.replace(' ','_')}_{plot_end_day_str.replace(' ','_')}_{st}_requested_{default_datum}.pkl"
        cache_success = False
        if os.path.exists(cache_filename):
            try:
                with open(cache_filename, 'rb') as f:  # Python 3: open(..., 'rb')
                    this_noaa_df, this_datum, station_data = pickle.load(f)
                    print(f'Existing obs data read from {cache_filename}')
                    cache_success = True
            except:
                print(f'Failed to read from Cache, regenerating cache')

        if not cache_success:
            try:
                if retrieve_method == 'noaa_coops':
                    station_data = noaa_coops.Station(st)
                elif retrieve_method == 'searvey': # searvey, convert into the same format as noaa_coops
                    station_data = COOPS_Station(int(st))
                    lon_lat = np.squeeze(np.array(station_data.location.coords))
                    setattr(station_data, 'lat_lon', {'lat': lon_lat[-1], 'lon': lon_lat[0]})
                else:
                    raise Exception(f'retrieve_method {retrieve_method} not supported')
            except:  # JSONDecodeError
                raise Exception("Got JSONDecodeError, possible unstable network from COOPS server")

            try:
                if retrieve_method == 'noaa_coops':
                    this_noaa_df = station_data.get_data(
                        begin_date=plot_start_day_str.replace("-", "")[:-9],
                        end_date=plot_end_day_str.replace("-", "")[:-9],
                        product="hourly_height", datum=default_datum, units="metric", time_zone="gmt"
                    )
                    this_noaa_df.reset_index(inplace=True)
                    this_noaa_df.columns = ['date_time', 'water_level', 'sigma', 'flags']
                    this_noaa_df = this_noaa_df.set_index('date_time')
                elif retrieve_method == 'searvey':
                    this_noaa_df = COOPS_Query(
                        station=st, start_date=plot_start_day_str, end_date=plot_end_day_str,
                        datum=default_datum, product='hourly_height',
                        units='metric', time_zone='gmt'
                    ).data
                    # convert into the same format as noaa_coops
                    this_noaa_df.reset_index(inplace=True)
                    this_noaa_df = this_noaa_df.rename(index='date_time')
                    substitute_col = {'t': 'date_time', 'v': 'water_level', 's': 'sigma', 'f': 'flags', 'q': 'QC'}
                    this_noaa_df = this_noaa_df.rename(columns=substitute_col)
                    this_noaa_df.set_index('date_time', inplace=True)
                else:
                    raise Exception(f'retrieve_method {retrieve_method} not supported')

                this_datum = default_datum
            except Exception:
                try:
                    this_noaa_df = station_data.get_data(
                        begin_date=plot_start_day_str.replace("-", "")[:-3],
                        end_date=plot_end_day_str.replace("-", "")[:-3],
                        product="hourly_height", datum="MSL", units="metric", time_zone="gmt"
                    )
                    this_noaa_df.reset_index(inplace=True)
                    this_noaa_df.columns = ['date_time', 'water_level', 'sigma', 'flags']
                    this_noaa_df = this_noaa_df.set_index('date_time')

                    this_datum = "MSL"
                except Exception:
                    this_noaa_df = None
                    this_datum = None

            with open(cache_filename, 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump([this_noaa_df, this_datum, station_data], f)

        noaa_df_list.append(this_noaa_df)
        datum_list.append(this_datum)
        station_data_list.append(station_data)
        st_info.append({'site_name': station_data.name, 'id': station_data.id, 'latitude': station_data.lat_lon['lat'], 'longitude': station_data.lat_lon['lon']})


    return [noaa_df_list, datum_list, st_info]

def get_usace_elev(stations=[], start_time_str='2022-04-11T00:00', end_time_str='2022-04-19T23:59'):
    import xml.etree.ElementTree as ET
    import requests
    import pytz
    '''
    Get gage height from USACE, using the following API:
    https://rivergages.mvr.usace.army.mil/watercontrol/webservices/rest/webserviceWaterML.cfc?method=RGWML&meth=getValues&location=01300&site=01300&variable=HG&beginDate=2022-04-11T00:00&endDate=2022-04-19T23:59&authToken=RiverGages
    '''
    cache_folder = os.path.realpath(os.path.expanduser('~/schism10/Cache/'))

    # check time format
    try:
        start_time = parse_date(start_time_str)[0]
        end_time = parse_date(end_time_str)[0]
    except:
        raise Exception(f'Invalid time format: {start_time_str}, {end_time_str}')
    # convert timezone from Louisiana local time to GMT
    gmt_tz = pytz.timezone('GMT')
    louisiana_tz = pytz.timezone('America/Chicago')
    start_time_str = gmt_tz.localize(start_time).astimezone(louisiana_tz).strftime('%Y-%m-%dT%H:%M')
    end_time_str = gmt_tz.localize(end_time).astimezone(louisiana_tz).strftime('%Y-%m-%dT%H:%M')


    df_list = []; st_info_list = []
    for i, station in enumerate(stations):
        # look for cache first
        print(f'Processing Station {i+1} of {len(stations)}: {station}')

        cache_filename = f"{cache_folder}/USACE_{start_time_str.replace(' ','_')}_{end_time_str.replace(' ','_')}_{station}.pkl"
        cache_success = False
        if os.path.exists(cache_filename):
            try:
                with open(cache_filename, 'rb') as f:  # Python 3: open(..., 'rb')
                    df, st_info = pickle.load(f)
                    print(f'Existing obs data read from {cache_filename}')
                    cache_success = True
            except:
                print(f'Failed to read from Cache, regenerating cache')

        if not cache_success:
            url = 'https://rivergages.mvr.usace.army.mil/watercontrol/webservices/rest/webserviceWaterML.cfc?' + \
                f'method=RGWML&meth=getValues&location={station}&site={station}' + \
                f'&variable=HG&beginDate={start_time_str}&endDate={end_time_str}&authToken=RiverGages'
            root = ET.fromstring(requests.get(url, verify=False).text)
            # namespace for waterML
            ns = {'wml': 'http://www.cuahsi.org/waterML/1.0/'}
            # Extract site information
            site_info = root.find('.//wml:sourceInfo', ns)
            site_name = site_info.find('wml:siteName', ns).text
            site_code = site_info.find('.//wml:siteCode', ns).text
            latitude = float(site_info.find('.//wml:latitude', ns).text)
            longitude = float(site_info.find('.//wml:longitude', ns).text)

            st_info = {'site_name': site_name, 'id': site_code, 'latitude': latitude, 'longitude': longitude}

            # Extract and print variable name
            variable_name = root.find('.//wml:variableName', ns).text
            print(f"Variable Name: {variable_name}")
            
            # Extract and iterate through all values
            values = root.findall('.//wml:value', ns)
            value_array = np.nan * np.ones(len(values), dtype=float)
            time_array = np.empty(len(values), dtype=object)
            for i, value in enumerate(values):
                time_array[i] = parse_date(value.get('dateTime'))[0]
                value_array[i] = float(value.text)

            # convert from feet to meter
            value_array = value_array * 0.3048
            # Create a pandas dataframe
            df = pd.DataFrame({'datetime': time_array, 'water_level': value_array})
            df.set_index('datetime', inplace=True)
            # convert to GMT
            df.index = df.index.tz_localize(louisiana_tz).tz_convert(gmt_tz)
            # set unreadable values to nan
            df.loc[df['water_level'] < -99, 'water_level'] = np.nan
            df.loc[df['water_level'] > 99, 'water_level'] = np.nan

            # save cache
            with open(cache_filename, 'wb') as f:
                pickle.dump([df, st_info], f)

        df_list.append(df)
        st_info_list.append(st_info)


    return df_list, st_info_list

def make_bp_from_station_json(
    bpfile='/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/station.bp',
    station_json='/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/station.json'
):

    _, st_info = get_timeseries_from_station_json(station_json=station_json)
    lon = [x['longitude'] for x in st_info]
    lat = [x['latitude'] for x in st_info]
    z = [0] * len(lon)
    st_id = [x['id'] for x in st_info]
    
    with open(bpfile, 'w') as f:
        f.write("\n")
        f.write(str(len(lon)) + "\n")
        for i, _ in enumerate(lon):
            f.write(f"{i+1} {lon[i]} {lat[i]} {z[i]} !{st_id[i]}\n")  # "!" is a separator for station id


def get_timeseries_from_station_json(case_name='Missi_ida', station_json='/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/station.json'):
    '''
    This is a sample function to get time series from different sources
    using the station json file as input.

    The outputs are obs and model time series in the
    format of pandas dataframe.
    '''

    # get station information from json
    with open(station_json) as d:
        dict = json.load(d)
    
    # get mixed sources
    stations = np.array(list(dict[case_name]['stations'].keys()))

    # coops
    coops_station_idx = [i for i, x in enumerate(dict[case_name]['stations'].values()) if x["source"] == 'COOPS']
    if coops_station_idx:
        [coops_obs, coops_datums, coops_st_info] = get_coops_elev(
            plot_start_day_str=dict[case_name]['plot_start_day_str'],
            plot_end_day_str=dict[case_name]['plot_end_day_str'],
            noaa_stations=stations[coops_station_idx],
            default_datum="NAVD"
        )
        # fix datums if the NAVD is not available
        for i, datum in enumerate(coops_datums):
            if datum != "NAVD":
                this_station = stations[coops_station_idx[i]]
                datum_shift = dict[case_name]['stations'][this_station]['to_NAVD88_feet']
                datum_shift *= 0.3048  # to meters
                coops_obs[i]['water_level'] += datum_shift
    
    # get usace
    usace_station_idx = [i for i, x in enumerate(dict[case_name]['stations'].values()) if x["source"] == 'USACE']
    if usace_station_idx:
        [usace_obs, usace_st_info] = get_usace_elev(
            stations=stations[usace_station_idx],
            start_time_str=dict[case_name]['plot_start_day_str'],
            end_time_str=dict[case_name]['plot_end_day_str']
        )
        # USCAE datums are local, need to be converted to NAVD88; USACE data is already converted from feet to meters
        for i, _ in enumerate(usace_obs):
            this_station = stations[usace_station_idx[i]]
            datum_shift = dict[case_name]['stations'][this_station]['to_NAVD88_feet']
            datum_shift *= 0.3048  # to meters
            usace_obs[i]['water_level'] += datum_shift
    
    # get usgs
    usgs_station_idx = [i for i, x in enumerate(dict[case_name]['stations'].values()) if x["source"] == 'USGS']
    if usgs_station_idx:
        usgs_obs, usgs_st_info = get_usgs_elev(
            station_ids=stations[usgs_station_idx],
            start_date=dict[case_name]['plot_start_day_str'],
            end_date=dict[case_name]['plot_end_day_str'],
        )
        # USGS datums may not be in NAVD88, need to be converted to NAVD88
        for i, _ in enumerate(usgs_obs):
            this_station = stations[usgs_station_idx[i]]
            datum_shift = dict[case_name]['stations'][this_station]['to_NAVD88_feet']
            datum_shift *= 0.3048
            usgs_obs[i]['water_level'] += datum_shift

    # assemble all obs data in order
    obs = []; st_info = []
    for i, station in enumerate(stations):
        if i in coops_station_idx:
            obs.append(coops_obs[coops_station_idx.index(i)])
            st_info.append(coops_st_info[coops_station_idx.index(i)])
        elif i in usace_station_idx:
            obs.append(usace_obs[usace_station_idx.index(i)])
            st_info.append(usace_st_info[usace_station_idx.index(i)])
        elif i in usgs_station_idx:
            obs.append(usgs_obs[usgs_station_idx.index(i)])
            st_info.append(usgs_st_info[usgs_station_idx.index(i)])


    # get model
    # mod = TimeHistory(
    #     file_name = dict[case_name]['elev_out_file'],
    #     start_time_str=dict[case_name]['model_start_day_str'],
    #     sec_per_time_unit=86400,
    #     columns=list(dict[case_name]['stations'].values())
    # )
    # mod.df.set_index('datetime', inplace=True)

    datums = ['NAVD' for _ in obs]

    return obs, st_info, datums


if __name__ == "__main__":
    make_bp_from_station_json()
    # plot_operation()
    # os.system("rm stats*.txt *png")
