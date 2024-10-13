import matplotlib.pyplot as plt
import noaa_coops
from searvey.coops import COOPS_Query, COOPS_Station
from schism_py_pre_post.Grid.Bpfile import Bpfile
from schism_py_pre_post.Shared_modules.test_modules import HYCOM
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Download.download_coops_elev import get_coops_water_level
from schism_py_pre_post.Utilities.util import vdatum_wrapper_pointwise, vdatum_preset
import pandas as pd
from datetime import timedelta
from datetime import datetime
from schism_py_pre_post.Shared_modules.obs_mod_comp import obs_mod_comp
from get_elev import get_obs_from_station_json
import numpy as np
from scipy.interpolate import griddata  # , interp2d
import math
import os
import glob
import pickle
import json
import copy
from pathlib import Path


def write_stat(stats, fname, readablity=2):
    stats_df = copy.deepcopy(stats)
    rearranged_cols = ['station_name', 'station_id', 'station_lon', 'station_lat', 'RMSE', 'MAE',
                       'Bias', 'CC', 'ubRMSE', 'Max_Obs', 'Max_Mod']
    # remove invalid obs before calculating mean
    stats_df = stats_df[stats_df['station_name'] != 'NA']
    stats_df = stats_df[stats_df['Max_Obs'] < 20]
    stats_df = stats_df[stats_df['Max_Mod'] < 20]

    # rearrange columns
    stats_df = stats_df[rearranged_cols]

    means = stats_df.iloc[:, 2:].mean(axis=0)
    mean_row = pd.DataFrame(columns=list(stats_df.columns),
                            data=(np.r_[np.nan, np.nan, means.values]).reshape(1, -1))
    stats_df = pd.concat([mean_row, stats_df], ignore_index=True)

    if readablity >= 1:
        stats_df.to_csv(
            Path(fname).with_suffix('.csv'), encoding='utf-8',
            index=False, na_rep='-9999', float_format='%.4f')

    if readablity >= 2:
        pd.options.display.float_format = '{:,.4f}'.format
        formatted_string = stats_df.to_string(index=False)
        with open(fname, 'w', encoding='utf-8') as f:
            f.write(formatted_string)

    # the first row (var names) and second row (values)
    # stats_df.iloc[:1, :].to_string(index=False).split('\n')

    return stats_df

def get_hycom_elev(noaa_stations=None, station_in_file=None, hycom_file=None):
    bp_df = Bpfile(station_in_file, cols=5).make_dataframe()
    if noaa_stations is not None:
        bp_df = bp_df[noaa_stations]

    st_lon = bp_df.loc['lon'].to_numpy()
    st_lat = bp_df.loc['lat'].to_numpy()

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


def get_hindcast_elev(
    model_start_day_str, noaa_stations=None, station_in_file=None,
    elev_out_file=None, sec_per_time_unit=1, station_in_subset=None
):
    '''
    Get model time series at noaa_stations,
    either from specified ids or from station_in_file
    '''
    if noaa_stations is not None:
        pass
    else:
        my_th = TimeHistory(elev_out_file, model_start_day_str, -9999, sec_per_time_unit=sec_per_time_unit)
        st_id = Bpfile(station_in_file, cols=5).make_dataframe().columns

        if station_in_subset is not None:
            my_th = my_th.export_subset(station_idx=station_in_subset)
            st_id = st_id[station_in_subset]
        
        model_df = my_th.df.set_index('datetime')
        model_df.columns = st_id

    return model_df

def datum_shift(original_df, datum_shift_file=None):
    # read datum_shift
    datum_shifts_df = pd.read_csv(datum_shift_file)
    datum_shifts_df['ID'] = datum_shifts_df['ID'].astype(str)

    # set shift on each station
    mod_shift = np.zeros((len(original_df.columns),), dtype=float)
    for i, station_id in enumerate(original_df.columns):
        if station_id in datum_shifts_df['ID'].values:
            mod_shift[i] = datum_shifts_df[datum_shifts_df['ID'] == station_id]['shift'].values[0]
    
    # expand shift to the same shape as mod
    mod_shift = np.tile(mod_shift, (original_df.shape[0], 1))
    # Apply shift. The shift is from NAVD to geoid, so we need to subtract it from mod (geoid to NAVD)
    shifted_df = original_df - mod_shift

    return shifted_df

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


def get_obs_elev(
    plot_start_day_str=None, plot_end_day_str=None, noaa_stations=None,
    default_datum='NAVD', cache_folder='/sciclone/schism10/feiye/Cache/',
    retrieve_method='noaa_coops',  # 'native', 'searvey' or 'noaa_coops'
):

    noaa_df_list = []
    datum_list = []
    station_data_list = []
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
            max_try = 5
            for ntry in range(max_try):
                try:
                    ntry += 1
                    if retrieve_method in ['native', 'noaa_coops']:
                        station_data = noaa_coops.Station(st)
                        break
                    elif retrieve_method == 'searvey': # searvey, convert into the same format as noaa_coops
                        station_data = COOPS_Station(int(st))
                        lon_lat = np.squeeze(np.array(station_data.location.coords))
                        setattr(station_data, 'lat_lon', {'lat': lon_lat[-1], 'lon': lon_lat[0]})
                        break
                    else:
                        raise Exception(f'retrieve_method {retrieve_method} not supported')
                except Exception as e:
                    print(f"Got Exception: {e} for {st} on try {ntry}, retrying")
                    if ntry >= max_try:
                        raise Exception("Got Exception: {e}, possible unstable network connecting to COOPS server")

            this_noaa_df = None; this_datum = None  # intialize
            for datum in [default_datum, 'MSL']:  # fall back to MSL, which is generally available
                print(f"Trying datum {datum}")
                try:
                    if retrieve_method == 'noaa_coops':
                        this_noaa_df = station_data.get_data(
                            begin_date=plot_start_day_str.replace("-", "")[:-9],
                            end_date=plot_end_day_str.replace("-", "")[:-9],
                            product="hourly_height", datum=datum, units="metric", time_zone="gmt"
                        )
                        this_noaa_df.reset_index(inplace=True)  # move t to a column
                        substitute_col = {'t': 'date_time', 'v': 'water_level', 's': 'sigma', 'f': 'flags', 'q': 'QC'}
                        this_noaa_df = this_noaa_df.rename(columns=substitute_col)
                        this_noaa_df.set_index('date_time', inplace=True)
                    elif retrieve_method == 'searvey':
                        this_noaa_df = COOPS_Query(
                            station=st, start_date=plot_start_day_str, end_date=plot_end_day_str,
                            datum=datum, product='hourly_height',
                            units='metric', time_zone='gmt'
                        ).data
                        # convert into the same format as noaa_coops
                        this_noaa_df.reset_index(inplace=True)
                        this_noaa_df = this_noaa_df.rename(index='date_time')
                        substitute_col = {'t': 'date_time', 'v': 'water_level', 's': 'sigma', 'f': 'flags', 'q': 'QC'}
                        this_noaa_df = this_noaa_df.rename(columns=substitute_col)
                        this_noaa_df.set_index('date_time', inplace=True)
                    elif retrieve_method == 'native':  # native method from this package
                        this_noaa_df = get_coops_water_level(
                            begin_date=plot_start_day_str.replace("-", "")[:-9],
                            end_date=plot_end_day_str.replace("-", "")[:-9],
                            datum=datum, station=st
                        )
                    else:
                        raise Exception(f'retrieve_method {retrieve_method} not supported')
                except Exception as e:
                    print(f"Failed to get data for {st} with datum {datum}: {e}")
                
                # clip data to the desired time range
                if this_noaa_df is not None:
                    this_noaa_df = this_noaa_df[plot_start_day_str:plot_end_day_str]
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


def plot_elev(
    obs_df_list, mod_df_all_stations, plot_start_day_str, plot_end_day_str,
              noaa_stations, datum_list, station_info,
              plot_name, iplot=True, subplots_shape=(10, None),
              fig_ax=None, line_styles=['r.', 'k'], font_size=6,
              nday_moving_average=0, demean=False,
              label_strs=[], label_rot=10):
    '''
    Plot time series and calculate stats.
    "obs_df_list" is a list of pd.DataFrame
    "mod_df" is a pd.DataFrame: 1st column is time; other columns are elevation, with the station name as the column name
    '''

    n_stations = len(obs_df_list)

    n_defined_shape = sum([x is not None for x in subplots_shape])
    if n_defined_shape == 0:
        n_subplot_row = 10
        n_subplot_col = math.ceil(n_stations / n_subplot_row)
    elif n_defined_shape == 2:
        n_subplot_row, n_subplot_col = subplots_shape
    else:
        n_subplot_row = math.ceil(n_stations / subplots_shape[1]) if subplots_shape[0] is None else subplots_shape[0]
        n_subplot_col = math.ceil(n_stations / subplots_shape[0]) if subplots_shape[1] is None else subplots_shape[1]

    if label_strs == []:
        label_strs = ['obs', 'model']

    # Plot
    if fig_ax is None:
        plt.rcParams.update({'font.size': font_size})
        fig, ax = plt.subplots(n_subplot_row, n_subplot_col, figsize=(25, 20))
        ax = ax.ravel()
    else:
        fig, ax = fig_ax

    n = 0
    stats = []; obs_maxs = []; mod_maxs = []
    for i, st in enumerate(noaa_stations):
        if st in mod_df_all_stations.columns:
            mod_df = mod_df_all_stations[st]
        else:
            raise Exception(f'cannot find {st} in mod_df_all_stations')
            mod_df = obs_df_list[i]['water_level']

        if obs_df_list[i] is None:
            obs_df = mod_df * np.nan
        else:
            obs_df = obs_df_list[i]['water_level']

        # clip data to the desired time range
        obs_df = obs_df.loc[(obs_df.index >= plot_start_day_str) & (obs_df.index <= plot_end_day_str)]
        mod_df = mod_df.loc[(mod_df.index >= plot_start_day_str) & (mod_df.index <= plot_end_day_str)]

        if (datum_list[i] != 'NAVD' or demean):
            obs_df = obs_df - obs_df.mean()
            mod_df = mod_df - mod_df.mean()

        my_comp = obs_mod_comp(obs=pd.DataFrame({'datetime': obs_df.index, 'value': obs_df}),
                               mod=pd.DataFrame({'datetime': mod_df.index, 'value': mod_df}))

        if nday_moving_average > 0:
            my_comp.get_moving_average(nday_avg=nday_moving_average)

        my_comp.cal_stats()

        if line_styles[0] is not None:
            ax[n].plot(my_comp.obs_df.time, my_comp.obs_df.value, line_styles[0], label=label_strs[0])
        if line_styles[1] is not None:
            ax[n].plot(my_comp.mod_df.time, my_comp.mod_df.value, line_styles[1], label=label_strs[1])

        obs_maxs.append(max(my_comp.obs_df['value']))
        mod_maxs.append(max(my_comp.mod_df['value']))
        stats.append(my_comp.stats_dict)

        datum_str = "NAVD 88" if datum_list[n] == "NAVD" else datum_list[n] 
        if obs_df_list[i] is None:
            title_str = f'{st}, no data\n'
        else:
            existing_title = ax[n].get_title()
            if existing_title == '':
                # title_str = f'{st}, {datum_str}, {station_info[n].name}\n'
                try:
                    title_str = f'{st}, {datum_str}; {station_info[n].name}\n; {my_comp.stats_str}'
                except AttributeError:
                    title_str = f"{st}, {datum_str}; {station_info[n]['site_name']}\n; {my_comp.stats_str}"
            else:
                title_str = f'{existing_title}\n; {my_comp.stats_str}'
        ax[n].title.set_text(title_str)

        ax[n].set_xlim([datetime.strptime(plot_start_day_str, "%Y-%m-%d %H:%M:%S"),
                        datetime.strptime(plot_end_day_str, "%Y-%m-%d %H:%M:%S")])
        ax[n].tick_params(labelrotation=label_rot)
        ax[n].set_ylim([np.nanmin(np.fmin(my_comp.obs_df.value, my_comp.mod_df.value)) - 0.5,
                        np.nanmax(np.fmax(my_comp.obs_df.value, my_comp.mod_df.value)) + 0.5])
        # ax[n].set_ylim([1, 9])
                       
        n = n + 1
    ax[0].legend()
    for i in range(n, n_subplot_col * n_subplot_row):
        ax[i].axis('off')

    fig.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=1.0)
    if plot_name is not None:
        plt.savefig(f'{plot_name}.png', dpi=400)
    if iplot:
        plt.show()

    # assemble station info
    try:
        station_name_list = [x.name for x in station_info]
    except AttributeError:
        station_name_list = [x['site_name'] for x in station_info]
    try:
        lon_list = [x.lat_lon['lon'] for x in station_info]
    except AttributeError:
        lon_list = [x['longitude'] for x in station_info]
    try:
        lat_list = [x.lat_lon['lat'] for x in station_info]
    except AttributeError:
        lat_list = [x['latitude'] for x in station_info]

    stat_df = pd.DataFrame({
        'station_id': noaa_stations,
        'station_name': station_name_list,
        'station_lon': lon_list,
        'station_lat': lat_list,
        'RMSE': [float(x['RMSE']) for x in stats],
        'ubRMSE': [float(x['ubRMSE']) for x in stats],
        'MAE': [float(x['MAE']) for x in stats],
        'Bias': [float(x['Bias']) for x in stats],
        'CC': [float(x['CC']) for x in stats],
        'Max_Obs': obs_maxs,
        'Max_Mod': mod_maxs
    })

    return [stat_df, [fig, ax]]


def plot_operation():
    os.system("rm stat*.txt *.png")
    # time
    plot_start_day_str = '2024-07-06 00:00:00'
    plot_end_day_str = '2024-07-10 23:00:00'

    forecast_end_day_str = '2024-07-08'  # corresponding to folder name
    fcst_folder = '/sciclone/schism10/feiye/STOFS3D-v7/Outputs/'  # '/sciclone/schism10/feiye/STOFS3D-v4/fcst/'
    remote_dir = fcst_folder  # sending saved plots to this folder

    # station presets>>
    station_bp_file = '/sciclone/schism10/feiye/STOFS3D-v6/fcst/run/station.in'
    noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id
    noaa_stations_groups = {
        'Florida': noaa_stations_all[:10],
        'Atlantic': noaa_stations_all[10:29],
        'GoME': noaa_stations_all[29:39],
        'GoMX_west': noaa_stations_all[41:60],
        'GoMX_east': noaa_stations_all[60:80],
        'Atlantic_inland1': noaa_stations_all[80:100],
        'Atlantic_inland2': noaa_stations_all[100:120],
        'GoMX_inland': noaa_stations_all[120:150],
        'Puerto_Rico': noaa_stations_all[150:164] + noaa_stations_all[39:41]
    }
    default_datums = {
        'Florida': 'NAVD',
        'Atlantic': 'NAVD',
        'GoME': 'NAVD',
        'GoMX_west': 'NAVD',
        'GoMX_east': 'NAVD',
        'Atlantic_inland1': 'NAVD',
        'Atlantic_inland2': 'NAVD',
        'GoMX_inland': 'NAVD',
        'Puerto_Rico': 'MSL'
    }

    # if mod is not in NAVD, apply datum shift
    datum_shift_file = '/sciclone/schism10/feiye/STOFS3D-v6/fcst/run/navd2xgeoid_shift.txt'


    # # test>
    # noaa_stations_groups = {
    #     'Test': ['9751364']
    # }
    # default_datum = {'Test': "NAVD"]

    # Ida>
    # noaa_stations_groups = {
    #     'Ida_4_stations': ['8760922', '8760721', '8761305', '8747437']
    # }
    # default_datums = {'Ida_4_stations': "NAVD"]

    # noaa_stations_groups = {
    #     'Ida_stations': ['8760922', '8760721', '8761724', '8761305', '8747437', '8761927', '8741533']
    # }
    # default_datums = {'Ida_stations': "NAVD"}

    # To compare with ADCIRC>
    # noaa_stations_groups = {
    #     'ADCIRC_GOMX': ["8779770", "8779748", "8775870", "8775792", "8775296", "8775283", "8775237",
    #                     "8774770", "8773701", "8773259", "8773037", "8771450", "8771341", "8771013",
    #                     "8770971", "8770822", "8770777", "8770613", "8770520", "8770475", "8764227",
    #                     "8764044", "8761927", "8761305", "8760922", "8760721", "8747437", "8737048",
    #                     "8735391", "8735180", "8729210", "8729108", "8728690", "8727520", "8726724",
    #                     "8726607", "8726520", "8726384", "8725520", "8725110"]
    # default_datum = {'ADCORC_GOMX': "NAVD"}

    stats = pd.DataFrame()
    for group_name in noaa_stations_groups:
        # mod = get_hycom_elev(
        #     noaa_stations=noaa_station_groups[group_name],
        #     station_in_file=station_bp_file,
        #     hycom_file='/sciclone/schism10/feiye/TEMP/GLBy0_expt_93.nc'
        # )

        # SCHISM's staout_1
        # mod = get_hindcast_elev(
        #     model_start_day_str='2021-05-01 00:00:00',
        #     noaa_stations=None,
        #     station_in_file="/sciclone/schism10/hyu05/NOAA_NWM/oper_3D/fcst/20210829/station.in",
        #     elev_out_file='/sciclone/schism10/feiye/ICOGS/RUN10a2/PostP/staout_1'
        # )

        # STOFS-3D's forecast
        mod = get_forecast_elev(
            plot_start_day_str=plot_start_day_str,
            forecast_end_day_str=forecast_end_day_str,  # may be different from plot_end_day_str
            fcst_folder=fcst_folder,
            station_in_file=station_bp_file,
            i_nowcast=False
        )

        # datum_shift for model
        if datum_shift_file != '' or datum_shift_file is not None:
            datum_shifts_df = pd.read_csv(datum_shift_file)
            datum_shifts_df['ID'] = datum_shifts_df['ID'].astype(str)

            # add shift
            mod_shift = np.zeros((len(mod.columns),), dtype=float)
            for i, station_id in enumerate(mod.columns):
                if station_id in datum_shifts_df['ID'].values:
                    mod_shift[i] = datum_shifts_df[datum_shifts_df['ID'] == station_id]['shift'].values[0]
            
            # expand shift to the same shape as mod
            mod_shift = np.tile(mod_shift, (mod.shape[0], 1))
            # Apply shift. The shift is from NAVD to geoid, so we need to subtract it from mod (geoid to NAVD)
            mod = mod - mod_shift

        # get obs
        [obs, datums, st_info] = get_obs_elev(
            plot_start_day_str=plot_start_day_str,
            plot_end_day_str=plot_end_day_str,
            noaa_stations=noaa_stations_groups[group_name],
            default_datum=default_datums[group_name]
        )

        # plot
        tmp = plot_elev(obs, mod, plot_start_day_str, plot_end_day_str, noaa_stations_groups[group_name], datums, st_info, group_name, iplot=False)
        stats = stats.append(tmp[0])

    stats.loc['mean'] = stats.iloc[:, 4:].copy().mean()
    stats.at['mean', 'station_id'] = 'all'
    stats.at['mean', 'station_name'] = 'all'
    stats.to_csv('stats_STOFS_3D.txt', encoding='utf-8', index=False)
    pass

    # upload to ccrm drive:
    os.system(f"scp 'stats_STOFS_3D.txt' *png {remote_dir}")

def plot_HYCOM():
    # time
    plot_start_day_str = '2021-10-02 00:00:00'
    plot_end_day_str = '2021-10-19 00:00:00'

    # station presets>>
    # stations, ICOGS v2 and v3>
    station_bp_file = '/sciclone/schism10/hyu05/NOAA_NWM/oper_3D/fcst/20210829/station.in'
    noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id
    noaa_stations_groups = {
        'Florida': noaa_stations_all[:10],
        'Atlantic': noaa_stations_all[10:29],
        'GoME': noaa_stations_all[29:39],
        'GoMX_west': noaa_stations_all[40:60],
        'GoMX_east': noaa_stations_all[60:80],
        'Atlantic_inland1': noaa_stations_all[80:100],
        'Atlantic_inland2': noaa_stations_all[100:120],
        'GoMX_inland': noaa_stations_all[120:150],
        'Puerto_Rico': noaa_stations_all[150:163]
    }
    default_datums = {
        'Florida': 'MSL',
        'Atlantic': 'MSL',
        'GoME': 'MSL',
        'GoMX_west': 'MSL',
        'GoMX_east': 'MSL',
        'Atlantic_inland1': 'MSL',
        'Atlantic_inland2': 'MSL',
        'GoMX_inland': 'MSL',
        'Puerto_Rico': 'MSL'
    }
    hycom_file = '/sciclone/schism10/feiye/TEMP/GLBy0_expt_93.nc'

    # # Pacific
    # plot_start_day_str = '2018-08-01 00:00:00'
    # plot_end_day_str = '2018-09-11 00:00:00'
    # station_bp_file = '/sciclone/schism10/feiye/Coastal_Act/MISC/station.in'
    # noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id
    # noaa_stations_groups = {
    #     'Hawaii': noaa_stations_all[:12],
    #     'West_Coast': noaa_stations_all[12:],
    # }
    # default_datums = {
    #     'Hawaii': 'MSL',
    #     'West_Coast': 'NAVD',
    # }
    # hycom_file = '/sciclone/schism10/feiye/Coastal_Act/MISC/SSH_1.nc'

    stats = pd.DataFrame()
    for group_name in noaa_stations_groups:
        mod = get_hycom_elev(
            noaa_stations=noaa_stations_groups[group_name],
            station_in_file=station_bp_file,
            hycom_file=hycom_file
        )

        # get obs
        [obs, datums, st_info] = get_obs_elev(
            plot_start_day_str=plot_start_day_str,
            plot_end_day_str=plot_end_day_str,
            noaa_stations=noaa_stations_groups[group_name],
            default_datum=default_datums[group_name]
        )

        # plot
        filename = f'{group_name}_{default_datums[group_name]}'
        tmp = plot_elev(obs, mod, plot_start_day_str, plot_end_day_str, noaa_stations_groups[group_name],
                        datums, st_info, plot_name=filename, iplot=False, nday_moving_average=0)
        stats = stats.append(tmp[0])

    stats.loc['mean'] = stats.iloc[:, -4:].copy().mean()
    stats.at['mean', 'station_id'] = 'all'
    stats.at['mean', 'station_name'] = 'all'
    stats.to_csv('stats.txt', encoding='utf-8', index=False)
    pass

    # upload to ccrm drive:
    os.system("scp 'stats.txt' *png $cdir/srv/www/htdocs/yinglong/feiye/ICOGS/Compare_with_HYCOM/")


def test():
    '''Test for combined stats.'''
    stats_csv_dict = {
        '03u': '/sciclone/home/feiye/s1/STOFS3D-v8/O03u/stats_LA_reforecast_repos_nontidal_dat.csv',
        '03z': '/sciclone/home/feiye/s1/STOFS3D-v8/O03z/stats_LA_reforecast_repos_nontidal_dat.csv',
        '04c': '/sciclone/home/feiye/s1/STOFS3D-v8/O04c/stats_LA_reforecast_repos_nontidal_dat.csv',
    }

    stats_dict = {}
    for key, val in stats_csv_dict.items():
        stats_dict[key] = pd.read_csv(val)

    stats_combined = pd.concat(list(stats_dict.values())).groupby(level=0).min()
    # recalculate mean
    stats_combined.iloc[0, 2:] = stats_combined.iloc[:, 2:].mean(axis=0)

    print('test done.')

    
def test_plot_usgs():
    '''Plot time series and calculate stats.'''

    case_name = 'LA_reforecast_repos_nontidal'
    station_json_fname = '/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/station.json'

    with open(station_json_fname, 'r', encoding='utf-8') as f:
        dict = json.load(f)

    output_dir = dict[case_name]['cdir']

    obs, st_info, datums = get_obs_from_station_json(case_name, station_json_fname)

    mod = TimeHistory(
        file_name=dict[case_name]['elev_out_file'],
        start_time_str=dict[case_name]['model_start_day_str'],
        sec_per_time_unit=86400,
        columns=list(dict[case_name]['stations'].keys())
    )
    mod.df.set_index('datetime', inplace=True)

    # shift for mod
    mod.df.iloc[:, :] += 0.0

    if dict[case_name]['default_datum'].lower() == 'navd88':
        pass
    elif dict[case_name]['default_datum'].lower() == 'xgeoid20b':
        # datum shift for mod, from xgeoid20b to NAVD
        st_lon = np.array([x['longitude'] for x in st_info])
        st_lat = np.array([x['latitude'] for x in st_info])
        st_shift = vdatum_wrapper_pointwise(
            x=st_lon, y=st_lat, z=np.zeros_like(st_lat),
            conversion_para=vdatum_preset['xgeoid20b_to_navd88'],
            print_info='\nConverting from xgeoid20b to NAVD88:\n'
        )
        for i in range(len(st_info)):
            mod.df.iloc[:, i] += st_shift[i]  # from xgeoid20b to NAVD
    elif dict[case_name]['default_datum'].lower() in ['msl', "de-mean"]:
        # demean both obs and mod, for all time
        # to demean for the plot period, use plot_elev() with demean=True
        mod.df.iloc[:, :] -= mod.df.mean(axis=0)
        for i, ob in enumerate(obs):
            if ob is not None:
                obs[i]['water_level'] -= obs[i]['water_level'].mean()
    else:
        raise ValueError(f"Unknown default_datum: {dict['default_datum']}")

    # plot
    # split stations into chunks, each with 25 stations
    stations_in_plots = [slice(i, i+25) for i in range(0, len(obs), 25)]
    total_stats_df = pd.DataFrame()
    for stations_in_plot in stations_in_plots:
        plot_name = (
            f"{output_dir}/ts_{case_name}_{dict[case_name]['elev_out_file'].split('.')[-1]}"
            f"_{stations_in_plot.start}_{stations_in_plot.stop}")
        stats_df, _ = plot_elev(
            obs_df_list=obs[stations_in_plot],
            mod_df_all_stations=mod.df.iloc[:, stations_in_plot],
            plot_start_day_str=dict[case_name]['plot_start_day_str'],
            plot_end_day_str=dict[case_name]['plot_end_day_str'],
            noaa_stations=list(dict[case_name]['stations'].keys())[stations_in_plot],
            datum_list=datums[stations_in_plot],
            demean=True,
            station_info=st_info[stations_in_plot],
            iplot=False, plot_name=plot_name,
            subplots_shape=(7, None), label_strs=['obs', 'model'],
            font_size=8,
        )
        total_stats_df = total_stats_df.append(stats_df)
        write_stat(stats_df, f"{plot_name}_stats.txt")

    write_stat(total_stats_df, f"{output_dir}/"
               f"stats_{case_name}_{dict[case_name]['elev_out_file'].split('.')[-1]}.txt")

    print('Main function done.')


if __name__ == "__main__":
    test_plot_usgs()
    plot_operation()
    test()
    print('Done.')
