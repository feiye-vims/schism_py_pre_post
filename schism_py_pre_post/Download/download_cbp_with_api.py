import pandas as pd
import numpy as np
import datetime


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


class GenericObsData():
    __slots__ = ['station_info', 'df']

    def __init__(self, station_info=None, df=None):
        self.station_info = station_info
        self.df = df


def GetCBP(stations=None, sample_time=None, varname='SALINITY'):
    var_dict = {
        'SALINITY': 83,
        'WTEMP': 123
    }
    if stations is None:
        stations = ['CB1.1', 'CB2.1', 'CB2.2', 'CB3.1', 'CB3.2', 'CB3.3C',
                    'CB4.1C', 'CB4.2C', 'CB4.3C', 'CB4.4', 'CB5.1', 'CB5.2',
                    'CB5.3', 'CB5.4', 'CB5.5', 'CB7.1S', 'CB7.2', 'CB7.3', 'CB7.4']
    if sample_time is None:
        sample_time = '2011-08-24'
    sample_time = pd.to_datetime(sample_time, format='%Y-%m-%d', errors='coerce')
    var_id = var_dict[varname]

    search_window = [sample_time-datetime.timedelta(days=360),
                     sample_time+datetime.timedelta(days=50)]
    search_window = [x.strftime(format='%m-%d-%Y') for x in search_window]

    # get station dictionary
    stations_info = pd.read_json('http://datahub.chesapeakebay.net/api.json/Station', dtype=str)
    station_ids = [stations_info.loc[stations_info['StationName'] == station]['StationId'].values[0]
                   for station in stations]
    url = 'http://datahub.chesapeakebay.net/api.JSON/WaterQuality/WaterQuality/' + \
          search_window[0] + '/' + search_window[1] + \
          '/0,1/2,4,6/12,13,15,35,36,2,3,7,33,34,23,24/' + \
          'Station/' + ','.join(station_ids) + f'/{var_id}'

    df_all = pd.read_json(url)
    df_all['SampleDate'] = pd.to_datetime(df_all['SampleDate'])
    cbp_dict = {}
    for station in stations:
        df = df_all[df_all['Station'] == station]
        lon = df['Longitude'].values[0]
        lat = df['Latitude'].values[0]

        dt = df['SampleDate'] - sample_time
        idx = np.where(abs(dt) == abs(dt).min())
        df = df.iloc[idx]

        df = df.sort_values('Depth')
        df = df[['Station', 'Depth', 'MeasureValue', 'SampleDate', 'Parameter']].copy()
        df = df.drop_duplicates(subset=["Depth"])

        cbp_dict[station] = GenericObsData(
            station_info={'station_name': station, 'lon': lon, 'lat': lat},
            df=df
        )

    return cbp_dict


def get_cbp_obs_for_stofs3d(outdir=None, sample_time='2015-09-18', varname=['sal']):
    '''
    Download from usgs using climata via api
    - List of param ids:
    - sample api url:
      http://datahub.chesapeakebay.net/api.JSON/WaterQuality/WaterQuality/
      8-18-2016/8-18-2021/0,1/2,4,6/12,13,15,35,36,2,3,7,33,34,23,24/Station/1150/83,123
    '''

    var_dict = {
        'sal': 'SALINITY',
        'tem': 'WTEMP',
    }

    for var in varname:
        cbp_var = var_dict[var]
        my_obs = GetCBP(sample_time=sample_time, varname=cbp_var)
        out = open(f'{outdir}/mean_{var}_xyz_{sample_time}', 'a')
        for station in my_obs:
            x = my_obs[station].station_info['lon']
            y = my_obs[station].station_info['lat']
            z = my_obs[station].df['MeasureValue'].values[-1]
            st = my_obs[station].station_info['station_name']
            out.write(f'{x} {y} {z} {st}\n')
        out.close()
        
    pass
