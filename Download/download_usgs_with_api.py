from climata.usgs import InstantValueIO, SiteIO  # pip install climata (mannually define a class "SiteIOIV" in __init__.py)
# from climata.usgs import DailyValueIO
# from climata.acis import StationDataIO
import pandas as pd
import numpy as np
import gsw
import os
# import subprocess
# from pandas.plotting import register_matplotlib_converters


class SiteIOIV(SiteIO):
    default_params = {
        'format': 'rdb,1.0',
        'outputDataTypeCd': 'iv',
        yield lst[i:i + n]


def get_usgs_stations_from_state(states=['VA', 'MD'], data_type='iv'):
    col_names = ['site_no', 'station_nm', 'state', 'dec_long_va', 'dec_lat_va', 'begin_date', 'end_date']
    col = [[], [], [], [], [], [], []]

    for state in states:
        site_params = SiteIOIV(state=state)
        for i, col_name in enumerate(col_names):
            if col_name == 'state':  # state is not in site_param
                col[i] += [state for site_param in site_params]
            else:
                col[i] += [eval(f"site_param.{col_name}") for site_param in site_params]
        print(f"{state} stations found: {len(site_params)}")

    df = pd.DataFrame(np.array(col).T, columns=col_names)
    print(f"Total stations found: {len(df['site_no'])}")
    return df

        'siteStatus': 'all',
    }


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):

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


def all_states_time_average(param_id=None, var=None, unit_conv=1, states=None,
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

    if states is None:
        states = ['ME', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'DE', 'PA', 'MD', 'VA', 'NC', 'SC', 'GA', 'FL', 'AL', 'MI', 'LA', 'TX']
    station_info_df = get_usgs_stations_from_state(states)
    station_ids = station_info_df["site_no"].to_list()
    stations_chunk_size = 200

    print(f'states: {states}\n')
    print(f'var: {var}\n')

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

    # write mean_val_xyz
    with open(outfilename, 'w') as fout:
        for data in total_data:
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


def download_single_station(
    station_id='04044755',
    param_id='00010', var='temperature', unit_conv=1,
    datelist=pd.date_range(start='2019-10-10', end='2019-11-10')
):
    data_chunk = InstantValueIO(start_date=datelist[0], end_date=datelist[-1],
                                station=station_id, parameter=param_id)
    data_list = Process_data_chunk(data_chunk)
    return data_list


if __name__ == "__main__":
    # input section
    vars = ['temperature', 'salinity', 'conductance']
    outdir = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/USGS_DATA/'
    start_date_str = '2012-10-15'
    end_date_str = '2012-10-16'
    # end input section

    # dict of param ids:
    usgs_var_dict = {
        "streamflow": {"id": "00060", "unit_conv": 1},
        "salinity": {"id": "00480", "unit_conv": 1},
        "gauge height": {"id": "00065", "unit_conv": 1},
        "temperature": {"id": '00010', "unit_conv": 1},
        "conductance": {"id": '00095', "unit_conv": 0.001},
    }

    outfilenames = []
    for var in vars:
        outfilenames.append(f'{outdir}/mean_{var}_xyz_{start_date_str}')
        all_states_time_average(
            var=var,
            param_id=usgs_var_dict[var]["id"],
            unit_conv=usgs_var_dict[var]["unit_conv"],
            datelist=pd.date_range(start=start_date_str, end=end_date_str),
            outfilename=outfilenames[-1]
        )

    cmd_str = f"cat {outfilenames[1]} {outfilenames[2]} > {outdir}/mean_salinity_cond_xyz_{start_date_str}"
    print(cmd_str); os.system(cmd_str)
    os.symlink(f"mean_salinity_cond_xyz_{start_date_str}", f"{outdir}/mean_sal_xyz_{start_date_str}")
    os.symlink(f"mean_temperature_cond_xyz_{start_date_str}", f"{outdir}/mean_tem_xyz_{start_date_str}")
