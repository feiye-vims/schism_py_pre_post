from functools import cache
from climata.usgs import InstantValueIO, SiteIO
import pandas as pd
import numpy as np
import gsw
import os
from datetime import datetime, timedelta
# import subprocess
# from pandas.plotting import register_matplotlib_converters
from pylib import schism_grid
import matplotlib.pyplot as plt
from schism_py_pre_post.Download.Data import ObsData, Station


# dict of param ids:
usgs_var_dict = {
    "streamflow": {"id": "00060", "unit_conv": 0.028316846592},
    "salinity": {"id": "00480", "unit_conv": 1},
    "gauge height": {"id": "00065", "unit_conv": 1},
    "temperature": {"id": '00010', "unit_conv": 1},
    "conductance": {"id": '00095', "unit_conv": 0.001},
}

ecgc_states = ['ME', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'DE', 'PA', 'MD', 'VA', 'NC', 'SC', 'GA', 'FL', 'AL', 'MI', 'LA', 'TX']


class SiteIOIV(SiteIO):
    default_params = {
        'format': 'rdb,1.0',
        'outputDataTypeCd': 'iv',
        'siteStatus': 'all',
    }


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def get_usgs_stations_from_state(states=['VA', 'MD'], data_type='iv', parameter=None):
    col_names = ['site_no', 'station_nm', 'state', 'dec_long_va', 'dec_lat_va', 'begin_date', 'end_date']
    col = [[], [], [], [], [], [], []]

    for state in states:
        site_params = SiteIOIV(state=state)
        for i, col_name in enumerate(col_names):
            if col_name == 'state':  # state is not in site_param
                if parameter is None:
                    new_col = [state for site_param in site_params]
                else:
                    new_col = [state for site_param in site_params if site_param.parm_cd==parameter]
                col[i] += new_col
                print(f"{state} stations found: {len(new_col)}")
            else:
                if parameter is None:
                    col[i] += [eval(f"site_param.{col_name}") for site_param in site_params]
                else:
                    col[i] += [eval(f"site_param.{col_name}") for site_param in site_params if site_param.parm_cd==parameter]

    df = pd.DataFrame(np.array(col).T, columns=col_names)
    print(f"Total stations found: {len(df['site_no'])}")
    return df


def Process_data_chunk(data_chunk):
    total_data = []
    for data in data_chunk:
        dates = [row[0] for row in data.data]
        values = [row[1] for row in data.data]
        df = pd.DataFrame({'date': dates, 'value': values})
        total_data.append(
            GenericObsData(
                station_info={'id': data.site_code,
                              'name': data.site_name,
                              'lon': data.longitude,
                              'lat': data.latitude,
                              'var_name': data.variable_name,
                              'var_code': data.variable_code,
                              'unit': data.unit},
                df=df
            )
        )
    return total_data


class GenericObsData():
    __slots__ = ['station_info', 'df']

    def __init__(self, station_info=None, df=None):
        self.station_info = station_info
        self.df = df

def download_stations(param_id=None, station_ids=[], datelist=pd.date_range(start='2000-01-01', end='2000-01-02')):
    stations_chunk_size = 200

    if type(station_ids) is np.ndarray:
        station_ids = station_ids.tolist()

    print(f'station ids: {station_ids}\n')

    # download
    total_data = []
    for i, station_ids_chunk in enumerate(chunks(station_ids, stations_chunk_size)):
        data_chunk = InstantValueIO(start_date=datelist[0], end_date=datelist[-1],
                                    station=station_ids_chunk,
                                    parameter=param_id)
        for data in data_chunk:
            dates = [row[0] for row in data.data]
            values = [row[1] for row in data.data]
            df = pd.DataFrame({'date': dates, 'value': values})
            total_data.append(
                GenericObsData(
                    station_info={'id': data.site_code,
                                  'name': data.site_name,
                                  'lon': data.longitude,
                                  'lat': data.latitude,
                                  'var_name': data.variable_name,
                                  'var_code': data.variable_code,
                                  'unit': data.unit},
                    df=df
                )
            )
        print(f'stations processed: {min(len(station_ids),(i+1)*stations_chunk_size)} of {len(station_ids)}')
    print(f'Number of stations with available data: {len(total_data)}')

    return total_data

def write_time_average(input_data, param_id=None, unit_conv=1, outfilename=None):
    # write mean_val_xyz
    with open(outfilename, 'w') as fout:
        for data in input_data:
            x = data.station_info['lon']
            y = data.station_info['lat']
            z = data.df['value'].mean()
            if param_id == '00095':  # convert conductance to salinity
                z = gsw.conversions.SP_from_C(z * unit_conv, 25, 0)
            if np.isnan(z):
                continue  # skip nan
            fout.write(f"{x} {y} {z} {data.station_info['id']}\n")

    # plt.plot(dates, flow)
    # plt.xlabel('Date')
    # plt.ylabel('Streamflow')
    # plt.title(series.site_name)
    # plt.xticks(rotation='vertical')
    # plt.show()
    # except:
    #   print(f'No salinity at station {sid}')

def all_states_time_average(param_id=None, unit_conv=1, states=None,
                            datelist=pd.date_range(start='2000-01-01', end='2000-01-02'),
                            outfilename=None):
    '''
    Download from usgs using climata via api
    - Specify a list of stations (by states) to speed up the downloading process
      The USGS api takes a maximum of about 650 stations (with shorter names like 01109403) at a time
    - The default states are the east coast states
    - List of param ids:
        streamflow = "00060"
        salinity = "00480", unit_conv = 1
        gauge height = "00065"
        temperature = '00010'
        conductance = '00095', unit_conv = 0.001
    '''
    
    if os.path.exists(outfilename):
        print(f'{outfilename} exists, skipping ...')
        return

    if states is None:
        states = ['ME', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'DE', 'PA', 'MD', 'VA', 'NC', 'SC', 'GA', 'FL', 'AL', 'MI', 'LA', 'TX']
    station_info_df = get_usgs_stations_from_state(states)
    station_ids = station_info_df["site_no"].to_list()

    if not os.path.exists(os.path.dirname(outfilename)):
        os.mkdir(os.path.dirname(outfilename))

    total_data = download_stations(param_id=param_id, station_ids=station_ids, datelist=datelist)
    write_time_average(input_data=total_data, param_id=param_id, unit_conv=unit_conv, outfilename=outfilename)


def download_single_station(
    station_id='04044755',
    param_id='00010',
    datelist=pd.date_range(start='2019-10-10', end='2019-11-10')
):
    data_chunk = InstantValueIO(start_date=datelist[0], end_date=datelist[-1],
                                station=station_id, parameter=param_id)
    data_list = Process_data_chunk(data_chunk)
    return data_list


def get_usgs_obs_for_stofs3d(outdir=None, start_date_str='2015-09-18', end_date_str=None, vars=None):
    if vars is None:
        vars = ['temperature', 'salinity', 'conductance']
    if end_date_str is None:
        end_date_str = (datetime.strptime(start_date_str, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
        
    if outdir is None:
        raise Exception('outdir not set')
    os.makedirs(outdir, exist_ok=True)

    outfilenames = []
    for var in vars:
        outfilenames.append(f'{outdir}/mean_{var}_xyz_{start_date_str}')
        all_states_time_average(
            param_id=usgs_var_dict[var]["id"],
            unit_conv=usgs_var_dict[var]["unit_conv"],
            datelist=pd.date_range(start=start_date_str, end=end_date_str),
            outfilename=outfilenames[-1]
        )

    cmd_str = f"cat {outfilenames[1]} {outfilenames[2]} > {outdir}/mean_salinity_cond_xyz_{start_date_str}"
    print(cmd_str); os.system(cmd_str)

    if os.path.exists(f"{outdir}mean_sal_xyz_{start_date_str}"):
        os.remove(f"{outdir}mean_sal_xyz_{start_date_str}")
    os.symlink(f"mean_salinity_cond_xyz_{start_date_str}", f"{outdir}/mean_sal_xyz_{start_date_str}")
    if os.path.exists(f"{outdir}/mean_tem_xyz_{start_date_str}"):
        os.remove(f"{outdir}/mean_tem_xyz_{start_date_str}")
    os.symlink(f"mean_temperature_xyz_{start_date_str}", f"{outdir}/mean_tem_xyz_{start_date_str}")

def convert_to_ObsData(total_data:list, cache_fname=''):
    '''
    Convert the downloaded "total_data" (see Sample 1 in __main__) to an old format used by some early scripts
    '''
    stations = []
    for data in total_data:
        my_station = Station(fname=None)
        my_station.id = data.station_info['id']
        my_station.lon = data.station_info['lon']
        my_station.lat = data.station_info['lat']
        my_station.df = data.df.set_index(pd.DatetimeIndex(pd.to_datetime(data.df["date"], utc=True)))
        stations.append(my_station)
    
    obs_data = ObsData(fname=None)
    obs_data.stations = stations
    obs_data.set_fID2id()
    obs_data.saved_file = cache_fname
    obs_data.save(obs_data.saved_file)

    return obs_data
    

if __name__ == "__main__":

    # Sample 1: Download discharge data for all stations in selected states
    # ------------------------- inputs --------------------------- 
    states = ['NJ', 'CT', 'NY']
    download_var = 'gauge height'  # check the dictionary at the beginning of this script
    outdir = '/sciclone/schism10/feiye/Test/'
    cache_fname = f"{outdir}/{'_'.join(states)}.pkl"
    start_date_str = '2023-01-01'
    end_date_str = '2023-01-21'
    # ------------------------- end inputs --------------------------- 
    
    # get station info
    station_info = get_usgs_stations_from_state(states=states, parameter=usgs_var_dict[download_var]["id"])
    station_lon = station_info['dec_long_va'].to_numpy()
    station_lat = station_info['dec_lat_va'].to_numpy()
    # get rid of invalid stations
    valid = (station_lon != '') & (station_lat != '')
    station_ids = station_info["site_no"].to_numpy()
    station_ids = station_ids[valid]
    station_lon = station_lon[valid].astype(float)
    station_lat = station_lat[valid].astype(float)

    # download data. 
    # "total_data" is the preferred format, which can be saved as *.pkl for later use
    # , or see Sample 3 for writing total_data to *.txt
    total_data = download_stations(
        param_id=usgs_var_dict[download_var]['id'],
        station_ids=station_ids,
        datelist=pd.date_range(start=start_date_str, end=end_date_str)
    )
    # note: the actual downloaded stations may be slightly different from the requested "station_ids" due to data availablility
    valid_stations = [data.station_info['id'] for data in total_data]

    # save data in an old obs format (ObsData) in *.pkl, because some earlier script reads that format
    convert_to_ObsData(total_data=total_data, cache_fname=cache_fname)
    pass

    '''
    # Sample 2: STOFS3D
    get_usgs_obs_for_stofs3d(vars=['gauge height'], outdir='/sciclone/schism10/feiye/STOFS3D-v4/Data/', start_date_str='2021-03-01', end_date_str='2021-03-20')

    # Sample 3: get stations inside hgrid and plot
    outdir = '/sciclone/schism10/feiye/STOFS3D-v4/Data/'

    station_info = get_usgs_stations_from_state(states=ecgc_states, parameter="00065")
    station_lon = station_info['dec_long_va'].to_numpy()
    station_lat = station_info['dec_lat_va'].to_numpy()

    valid = (station_lon != '') & (station_lat != '')
    station_ids = station_info["site_no"].to_numpy()
    station_ids = station_ids[valid]
    station_lon = station_lon[valid].astype(float)
    station_lat = station_lat[valid].astype(float)

    gd = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/hgrid.ll.pkl')  # gd.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/hgrid.ll.pkl')
    in_grid_idx = gd.inside_grid(pxy=np.c_[station_lon, station_lat])
    in_grid_idx = in_grid_idx.astype(bool)
    in_grid_station_ids = station_ids[in_grid_idx].tolist()

    total_data = download_stations(
        param_id=usgs_var_dict['gauge height']['id'],
        station_ids=in_grid_station_ids,
        datelist=pd.date_range(start='2021-03-01', end='2021-03-21')
    )
    valid_stations = [data.station_info['id'] for data in total_data]

    station_info = station_info.drop_duplicates(subset='site_no')
    station_info.iloc[np.where(station_info.site_no.isin(valid_stations))].to_csv(f'{outdir}/station_info.txt')

    nstations = len(total_data)
    plt.rcParams.update({'font.size': 10})
    valid_stations = []
    for i, total_data_chunk in enumerate(chunks(total_data, 25)):
        fig, ax = plt.subplots(5, 5, figsize=(25, 20))
        ax = ax.ravel()
        for n, data in enumerate(total_data_chunk):
            ax[n].plot(data.df['date'], data.df['value'])
            ax[n].title.set_text(data.station_info['id'])
            data.df.to_csv(f"{outdir}/{data.station_info['id']}.txt", index=False)
        plt.savefig(f'{outdir}/Chunk_{i}.png')

    # data = download_single_station(station_id='04044755', param_id='00060',
    #                                unit_conv='0.028316846592',
    #                                datelist=pd.date_range(start='2022-01-01', end='2022-02-10'))
'''
    pass