import matplotlib.pyplot as plt
import noaa_coops as noaa
from schism_py_pre_post.Grid.Bpfile import Bpfile
from schism_py_pre_post.Shared_modules.test_modules import HYCOM
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
import pandas as pd
from datetime import timedelta
from datetime import datetime
from schism_py_pre_post.Shared_modules.obs_mod_comp import obs_mod_comp
import numpy as np
from scipy.interpolate import griddata  # , interp2d
import math
import os
import glob
import pickle


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


def get_hindcast_elev(model_start_day_str, noaa_stations=None, station_in_file=None, elev_out_file=None, sec_per_time_unit=1):
    '''
    if noaa_stations = None, then read stations from station.in
    '''
    st_id = Bpfile(station_in_file, cols=5).make_dataframe().columns

    my_th = TimeHistory(elev_out_file, model_start_day_str, -9999, sec_per_time_unit=sec_per_time_unit)
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


def get_obs_elev(plot_start_day_str, plot_end_day_str, noaa_stations, default_datum='NAVD', cache_folder='./'):

    cache_filename = f"{cache_folder}/COOPS_{plot_start_day_str.replace(' ','_')}_{plot_end_day_str.replace(' ','_')}_{'.'.join(noaa_stations)}_{default_datum}.pkl"

    if os.path.exists(cache_filename):
        with open(cache_filename, 'rb') as f:  # Python 3: open(..., 'rb')
            noaa_df_list, datum_list, station_data = pickle.load(f)
            print(f'Existing obs data read from {cache_filename}')
    else:  # download NOAA COOPS obs
        noaa_df_list = []
        datum_list = []
        station_data = []
        for i, st in enumerate(noaa_stations):
            print(f'Processing Station {i+1} of {len(noaa_stations)}: {st}')
            noaa_data = noaa.Station(st)
            station_data.append(noaa_data)
            try:
                noaa_df_list.append(noaa_data.get_data(begin_date=plot_start_day_str.replace("-", "")[:-3],
                                    end_date=plot_end_day_str.replace("-", "")[:-3],
                                    product="water_level", datum=default_datum, units="metric", time_zone="gmt", interval='h'))
                datum_list.append(default_datum)
            except Exception:
                try:
                    noaa_df_list.append(noaa_data.get_data(begin_date=plot_start_day_str.replace("-", "")[:-3],
                                        end_date=plot_end_day_str.replace("-", "")[:-3],
                                        product="water_level", datum="MSL", units="metric", time_zone="gmt", interval='h'))
                    datum_list.append("MSL")
                except Exception:
                    noaa_df_list.append(None)
                    datum_list.append(None)
        with open(cache_filename, 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([noaa_df_list, datum_list, station_data], f)

    return [noaa_df_list, datum_list, station_data]


def plot_elev(obs_df_list, mod_df_all_stations, plot_start_day_str, plot_end_day_str,
              noaa_stations, datum_list, station_info, plot_name, iplot=True, subplots_shape=(10, None),
              fig_ax=None, line_styles=['r.', 'k'], shift=0, nday_moving_average=0, label_strs=[], label_rot=10):
    '''
    Plot time series and calculate stats.
    "obs_df_list" is a list of pd.DataFrame
    "mod_df" is a pd.DataFrame: 1st column is time; other columns are elevation, with the station name as the column name
    '''
    font_size = 10

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
        title_str = None

        if st in mod_df_all_stations.columns:
            mod_df = mod_df_all_stations[st]
            mod_df.values[:] += shift
        else:
            mod_df = obs_df_list[i]['water_level']

        if obs_df_list[i] is None:
            title_str = f'{i+1}: {st}, no data\n'
            obs_df = mod_df * np.nan
        else:
            obs_df = obs_df_list[i]['water_level']

        mask = (mod_df.index > plot_start_day_str) & (mod_df.index <= plot_end_day_str)
        mod_df = mod_df.loc[mask]

        # if False:
        if (datum_list[i] == 'MSL'):
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

        if title_str is None:
            existing_title = ax[n].get_title()
            datum_str = "NAVD 88" if datum_list[n] == "NAVD" else datum_list[n] 
            if existing_title == '':
                # title_str = f'{st}, {datum_str}, {station_info[n].name}\n'
                title_str = f'{st}, {datum_str}; {station_info[n].name}\n; {my_comp.stats_str}'
            else:
                title_str = f'{existing_title}; {my_comp.stats_str}'
            ax[n].title.set_text(title_str)

        ax[n].set_xlim([datetime.strptime(plot_start_day_str, "%Y-%m-%d %H:%M:%S"),
                        datetime.strptime(plot_end_day_str, "%Y-%m-%d %H:%M:%S")])
        ax[n].tick_params(labelrotation=label_rot)
        ax[n].set_ylim([np.nanmin(np.fmin(my_comp.obs_df.value, my_comp.mod_interp_df.value)) - 0.5,
                        np.nanmax(np.fmax(my_comp.obs_df.value, my_comp.mod_interp_df.value)) + 0.5])
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

    stat_df = pd.DataFrame({
        'station_id': noaa_stations,
        'station_name': [x.name for x in station_info],
        'station_lon': [x.lat_lon['lon'] for x in station_info],
        'station_lat': [x.lat_lon['lat'] for x in station_info],
        'RMSE': [float(x['RMSE']) for x in stats],
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
    plot_start_day_str = '2022-05-02 00:00:00'
    plot_end_day_str = '2022-05-13 00:00:00'

    forecast_end_day_str = '2022-05-11 00:00:00'  # corresponding to folder name
    fcst_folder = '/sciclone/schism10/hyu05/NOAA_NWM/oper_3D/fcst/'
    remote_dir = '$cdir/srv/www/htdocs/yinglong/feiye/ICOGS/STOFS-3D_fcst/2022_05_10/'

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
        'Florida': 'NAVD',
        'Atlantic': 'NAVD',
        'GoME': 'NAVD',
        'GoMX_west': 'NAVD',
        'GoMX_east': 'NAVD',
        'Atlantic_inland1': 'NAVD',
        'Atlantic_inland2': 'MSL',
        'GoMX_inland': 'NAVD',
        'Puerto_Rico': 'MSL'
    }

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

    stats.loc['mean'] = stats.iloc[:, 4:].mean()
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

    stats.loc['mean'] = stats.iloc[:, -4:].mean()
    stats.at['mean', 'station_id'] = 'all'
    stats.at['mean', 'station_name'] = 'all'
    stats.to_csv('stats.txt', encoding='utf-8', index=False)
    pass

    # upload to ccrm drive:
    os.system("scp 'stats.txt' *png $cdir/srv/www/htdocs/yinglong/feiye/ICOGS/Compare_with_HYCOM/")


if __name__ == "__main__":
    plot_operation()
    os.system("rm stats*.txt *png")
