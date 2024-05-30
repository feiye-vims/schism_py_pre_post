import pandas as pd
import time

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
    for i in range(3):
        print(f"Attempt {i+1}")
        response = requests.get(api_url, params=params)

        if response is not None:
            data = response.json()

            if 'error' in data.keys():
                print(f"Error: {data['error']['message']}")
                break

            if 'data' not in data.keys() or len(data['data']) == 0:
                print(f"No data available for the specified parameters, keep trying.")
                sleep(1)
                continue
            else:
                return data  # found data, return

        else:  # retry if status code is not 200
            time.sleep(1)  # Sleep for 1 second to avoid hitting the API too frequently

    return None

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

if __name__ == "__main__":
    # Example usage
    df = get_coops_water_level(begin_date='20240312', end_date='20240326', station='8548989', datum='NAVD', time_zone='gmt', units='metric')
    print(df.head())

    # Save data to CSV
    df.to_csv('water_level_data.csv')
