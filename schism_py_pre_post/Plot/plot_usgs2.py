"""
Model-data comparison for USGS gages
"""


import numpy as np
import pandas as pd
import requests
import json
from tqdm import tqdm


def get_usgs_station_info(site_id):
    '''
    Get USGS site information for a given site ID
    '''

    # Define the correct API endpoint
    base_url = "https://waterservices.usgs.gov/nwis/site/"

    # Set the parameters (note no JSON format is available here; RDB is used)
    params = {
        "format": "rdb",  # RDB is a supported format
        "sites": site_id  # USGS site number
    }

    try:
        # Send the API request
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            print(f"Error retrieving data: {response.text}")
            return None

        # Parse the RDB response (lines beginning with '#' are comments)
        lines = response.text.splitlines()
        data_lines = [line for line in lines if not line.startswith("#")]

        # Extract the header and data
        if len(data_lines) < 2:
            print("No site data available for the provided site ID.")
            return
        header = data_lines[0].split("\t")
        data = data_lines[2].split("\t")

        # Combine header and data into a dictionary for readability
        site_info = dict(zip(header, data))

    except Exception as e:
        print(f"Error retrieving data: {e}")
    
    return site_info


if __name__ == "__main__":
    # Example usage
    case_name = 'LA_reforecast_repos_nontidal'  # 'v8'
    station_json_fname = ('/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/'
                          'station.json')

    with open(station_json_fname, 'r', encoding='utf-8') as f:
        plot_dict = json.load(f)[case_name]
    station_dict = plot_dict['stations']

    # get station info from USGS api
    stations = np.array(list(station_dict.keys()))
    station_needing_datum = []
    for station in tqdm(stations):
        station_dict[station]['to_NAVD88_feet'] = None  # initialize

        site_info = get_usgs_station_info(station)

        if site_info is None:
            station_dict[station]['to_NAVD88_feet'] = None
            station_needing_datum.append(station)
            print(f"Failed to retrieve station info for {station}.")
            continue

        if site_info['alt_datum_cd'] == 'NAVD88':
            station_dict[station]['to_NAVD88_feet'] = float(site_info['alt_va'])
        else:
            station_needing_datum.append(station)
            print(f"Failed to retrieve station datum for {station}.")

    # read additional corrections from csv
    station_datum_correction = pd.read_csv(
        '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/usgsCorrectionsAtl1.csv'
    )
    # strip "US" from Gauge ID
    station_datum_correction['Gauge ID'] = station_datum_correction['Gauge ID'].str[2:]
    for station in station_needing_datum.copy():
        correction = station_datum_correction[station_datum_correction['Gauge ID'] == station]
        if correction.empty:
            print(f"No correction found for station {station}.")
            continue
        if len(correction) > 1:
            print(f"Multiple corrections found for station {station}, skipping.")
            continue

        datum = correction['Datum'].values[0]
        correction1 = correction['Correction(m)'].values[0]
        correction2 = correction['NGVD29 to NAVD88 correction (m) ADD to USGS data'].values[0]
        if correction1 in (' ', '-') or correction2 in (' ', '-') or datum == 'Missing Datum':
            print(f"No correction value found for station {station}.")
            continue

        station_dict[station]['to_NAVD88_feet'] = ((
            float(correction['Correction(m)'].values[0])
            + float(correction['NGVD29 to NAVD88 correction (m) ADD to USGS data'].values[0])
        ) * 3.28084)
        station_needing_datum.remove(station)

    # save new station json
    plot_dict['stations'] = station_dict
    plot_dict['default_datum'] = 'NAVD88'

    new_station_json_fname = station_json_fname.replace('.json', '3.json')
    with open(new_station_json_fname, 'w', encoding='utf-8') as f:
        json.dump(plot_dict, f, indent=4)
    
    print('Done!')

    # obs, st_info, datums = get_obs_from_station_json(case_name, station_json_fname)