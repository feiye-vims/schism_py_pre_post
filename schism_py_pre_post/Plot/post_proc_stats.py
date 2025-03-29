'''Post process the statistics of the simulation'''


import numpy as np
import pandas as pd


stat_files = {
    'R11a': '/sciclone/home/feiye/s1/STOFS3D-v8/O11a/stats_v8_dat.csv',
    'R11b': '/sciclone/home/feiye/s1/STOFS3D-v8/O11b/stats_v8_dat.csv',
    'RUN16_v6': '/sciclone/home/feiye/s1/STOFS3D-v8/O16_v6/stats_v8_2018_dat.csv',
    # 'R13_v7': '/sciclone/home/feiye/s1/STOFS3D-v8/O13_v7/USGS_2018_July/stats_v8_2018_dat.csv',
    'R20b': '/sciclone/home/feiye/s1/STOFS3D-v8/O20b/stats_v8_2018_dat.csv',
}

neglected_stations = [
    2300300,  # South Fork Little Manatee River at Wimauma
    2309220,
    2309425,
    2309447,
    2313230,
    2313231,
    2330150,  # outside domain
    253044080555900,
    253047080555600,
    264514080550700,
    7375222,
    7380401, 7380500, 7381000,
]
stats = []
for run_id, stat_file in stat_files.items():
    print(stat_file)
    stat = pd.read_csv(stat_file)

    # filter out suspicious records
    invalid = np.zeros(stat.shape[0], dtype=bool)
    # invalid[:75] = True
    # for station_id in neglected_stations:
    #     invalid = invalid | (stat['station_id'] == station_id)
    invalid = invalid | (stat['station_id'] == -9999)
    # invalid = invalid | (stat['RMSE'] > 10)
    invalid = invalid | (np.abs(stat['Bias']) > 5)
    invalid = invalid | (stat['Max_Obs'] > 20)
    invalid = invalid | (stat['Max_Mod'] > 20)
    # invalid = invalid | (np.abs(stat['ubRMSE']) > 5)

    print(f'Invalid records: {invalid.sum()} of {stat.shape[0]}')
    stat = stat[~invalid]
    stats.append(stat)


# find common station ids
common_stations = set(stats[0]['station_id'])
for stat in stats[1:]:
    common_stations = common_stations.intersection(stat['station_id'])

# calculate mean stats from common stations
for run_id, stat in zip(stat_files.keys(), stats):
    stat = stat[stat['station_id'].isin(common_stations)]
    print(f'----- Run: {run_id }----- ')
    print(stat[:].mean(numeric_only=True))

print('Done')