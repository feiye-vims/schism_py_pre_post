import pandas as pd
import json

def read_json_via_api(**params):
    '''
    download json using api:

    Example url:
    https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=20240301&end_date=20240301&station=8725520&product=water_level&datum=NAVD&time_zone=gmt&units=metric&application=DataAPI_Sample&format=json

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
    import requests
    api_url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

    # Make a GET request to the API
    # try 10 times
    data = None
    for i in range(10):
        response = requests.get(api_url, params=params)
        if response.status_code == 200:
            data = response.json()
            break

    if data is None:
        print("Failed to retrieve data. Status code:", response.status_code)
    
    return data

def get_coops_station_metadata(station='8725520'):
    pass
        

def get_coops_water_level(begin_date='20240301', end_date='20240301', station='8725520', datum='NAVD', time_zone='gmt', units='metric'):
    '''
    Get water level data from NOAA CO-OPS API at a specific station

    '''
    data = read_json_via_api(
        begin_date=begin_date, end_date=end_date, station=station,
        product="water_level", datum=datum, time_zone=time_zone, units=units,
        application="DataAPI_Sample", format="json"
    )

    # Extract data points
    data_points = data['data']

    # Create lists to store datetime and water level values
    datetimes = []
    water_levels = []

    # Iterate through data points and extract datetime and water level values
    for point in data_points:
        datetimes.append(point['t'])
        water_levels.append(float(point['v']))  # Convert water level to float

    # Create DataFrame
    df = pd.DataFrame({'date_time': datetimes, 'water_level': water_levels})
    df['date_time'] = pd.to_datetime(df['date_time'])
    df.set_index('date_time', inplace=True)

    return df
