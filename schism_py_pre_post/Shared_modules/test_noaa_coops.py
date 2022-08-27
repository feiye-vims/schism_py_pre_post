# %%
import noaa_coops as nc
import numpy as np
from pandas import DataFrame, read_csv
import copy


# %% find datum of specified stations
# FEET_TO_METERS = 0.3048
requested_datums = ['STND', 'MSL', 'NAVD88', 'MLLW']
other_fields = ['units', 'accepted']

# get st ids from csv
sta_file = '/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/AWIPS/SCHISM_SAMPLES/SHEF/stations_test.csv'
out_file = '/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/AWIPS/SCHISM_SAMPLES/SHEF/stations_test_datum_info.csv'
df0 = read_csv(sta_file)
st_ids = df0['station_id']

# get st ids from station.in
sta_file = '/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/COOPS_id2name/station.in'
out_file = '/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/COOPS_id2name/station_id_name.txt'
with open(sta_file) as f:
    f.readline()
    f.readline()
    st_ids = [line.strip().split('!')[1] for line in f]

# st_names = [nc.Station(id).name for id in st_ids]
# st_names = [st_name.replace('\"','') for st_name in st_names]

# %%
st_info = {}
for requested_datum in requested_datums:
    st_info[requested_datum] = np.nan * np.ones(len(st_ids))
for field in other_fields:
    st_info[field] = [None] * len(st_ids)
st_info['name'] = [None] * len(st_ids)
st_info['lon'] = np.nan * np.ones(len(st_ids))
st_info['lat'] = np.nan * np.ones(len(st_ids))
st_info['units'] = np.repeat(['N/A      '], len(st_ids))  # more spaces to allow longer strings

for i, st_id in enumerate(st_ids):
    try:
        st = nc.Station(st_id)
    except KeyError:
        print(f'warning: station {i+1}: {st_id} not found')

    print(f"station {i+1} of {len(st_ids)}: {st_id}")

    if hasattr(st, 'datums'):
        if st.datums['datums'] is None:
            continue
        for field in other_fields:
            st_info[field][i] = st.datums[field]

        for datum in st.datums['datums']:
            if datum['name'] in requested_datums:
                st_info[datum['name']][i] = datum['value']

    st_info['lon'][i], st_info['lat'][i], st_info['units'][i], st_info['name'][i] =\
    st.lat_lon['lon'], st.lat_lon['lat'], st.datums['units'], st.name

# make dataframe
df = DataFrame(
    columns=['station_id', 'station_name', 'lon', 'lat'] + requested_datums + ['units', 'accepted']
)
df['station_id'] = st_ids
df['station_name'] = st_info['name']
df['lon'], df['lat'] = st_info['lon'], st_info['lat']
for requested_datum in requested_datums:
    df[requested_datum] = st_info[requested_datum]
df['units'] = st_info['units']
df['accepted'] = st_info['accepted']
df.to_csv(out_file, index=False, na_rep='NaN', sep=';')

pass
