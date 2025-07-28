"""
Harmonic analysis of a time series and plot model-data comparison
"""


import os
import pickle
import numpy as np
import pandas as pd
import json
from tqdm import tqdm
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import geopandas as gpd

from utide import solve
from plot_elev import get_hindcast_elev, get_obs_elev
from schism_py_pre_post.Shared_modules.obs_mod_comp import obs_mod_comp
from pylib import read_schism_bpfile

cache_folder = os.path.realpath(os.path.expanduser('~/schism10/Cache/'))
common_tidal_constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']


def get_stations_from_groups(grouping_polygon_shpfile=None, station_bp_file=None):
    stations_groups = {}
    station_bp = read_schism_bpfile(station_bp_file)

    station_points_gdf = gpd.GeoDataFrame(
        {'Station_Id': station_bp.station},
        geometry=gpd.points_from_xy(station_bp.x, station_bp.y)
    )

    group_polys = gpd.read_file(grouping_polygon_shpfile)

    for i, group in enumerate(group_polys['Group_Name']):
        print('Group: ', group)
        # find stations in this polygon
        stations_in_group = station_points_gdf[station_points_gdf.geometry.within(group_polys.geometry[i])]
        stations_groups[group] = stations_in_group['Station_Id'].tolist()
        
    return stations_groups


def plot_ha(event=None, color_groups=None, cache: bool = False):
    '''
    Compare harmonic analysis of model and obs for each station.

    Parameters
    - event: dict, event information, which is read from a json file
    - color_groups: list of colors for each group (e.g., [Gulf_Coast, East_Coast],
      ['yellow', 'blue'])
    - cache: bool, whether to read from cache or not
      The cache file saves the harmonic analysis results for each station.
    '''

    cache_file = f'{event["cdir"]}/{event["plot_start_day_str"]}_{event["plot_end_day_str"]}_ha.pkl'

    if not cache:
        mod = get_hindcast_elev(
            model_start_day_str=event['model_start_day_str'],
            noaa_stations=None,
            station_in_file=event['station_bp_file'],
            elev_out_file=event['elev_out_file'],
            station_in_subset=np.arange(0, 164),  # the first 164 stations, which are COOPS stations
        )
        stations = np.array(mod.columns)

        [obs, datums, st_info] = get_obs_elev(
            plot_start_day_str=event['plot_start_day_str'],
            plot_end_day_str=event['plot_end_day_str'],
            noaa_stations=stations,
            cache_folder=cache_folder,
            retrieve_method='noaa_coops'
        )
        valid_idx = np.full((len(obs), ), True, dtype=bool)
        for i, ob in enumerate(obs):
            if ob is None:
                valid_idx[i] = False
                continue

            dt = ob.index[1] - ob.index[0]
            if len(ob) != 1 + (
                pd.Timestamp(event['plot_end_day_str']) - pd.Timestamp(event['plot_start_day_str'])
            ).total_seconds() / dt.total_seconds():
                valid_idx[i] = False

        # set color for each group
        color_idx = np.zeros(len(stations))
        if color_groups is not None:
            for i, [group_name, group_stations] in enumerate(color_groups.items()):
                for j, st in enumerate(stations):
                    if st in group_stations:
                        color_idx[j] = i


        # harmonic analysis
        ha_obs_results = np.empty(len(stations), dtype=object)
        ha_mod_results = np.empty(len(stations), dtype=object)
        for i, (st, info) in enumerate(tqdm(zip(stations, st_info), total=len(stations))):
            if not valid_idx[i]:
                print('Missing obs for station: ', st)
            else:
                mod_df = pd.DataFrame({'datetime': mod.index, 'value': np.squeeze(mod[st]).values})
                obs_df = pd.DataFrame({'datetime': obs[i].index, 'value': obs[i]['water_level']})
                my_comp = obs_mod_comp(obs_df, mod_df)
                ha_obs_results[i] = solve(
                    t=my_comp.obs_df.index, u=my_comp.obs_df['value'],
                    lat=info.lat_lon['lat'], constit=common_tidal_constituents)
                ha_mod_results[i] = solve(
                    t=my_comp.obs_df.index, u=my_comp.mod_interp_df['value'],
                    lat=info.lat_lon['lat'], constit=common_tidal_constituents)

        # Compute complex error for 8 common tidal constituents at each station
        obs_complex = np.full((len(stations), len(common_tidal_constituents)), np.nan, dtype=complex)
        mod_complex = np.full((len(stations), len(common_tidal_constituents)), np.nan, dtype=complex)
        complex_error = np.full((len(stations), len(common_tidal_constituents)), np.nan, dtype=float)
        complex_error_all_const = np.full((len(stations), ), np.nan, dtype=float)
        for i, (st, info, ob, md) in enumerate(zip(stations, st_info, ha_obs_results, ha_mod_results)):
            if not valid_idx[i]:
                continue
            for j, cons in enumerate(common_tidal_constituents):
                ob_idx = np.argwhere(ob.name == cons)[0][0]
                md_idx = np.argwhere(md.name == cons)[0][0]

                obs_complex[i, j] = ob.A[ob_idx] * np.exp(1j * np.deg2rad(ob.g[ob_idx]))
                mod_complex[i, j] = md.A[md_idx] * np.exp(1j * np.deg2rad(md.g[md_idx]))
                complex_error[i, j] = 1 / np.sqrt(2) * np.abs(
                    obs_complex[i, j] - mod_complex[i, j])
            # all constituents
            complex_error_all_const[i] = np.sqrt(np.sum(complex_error[i] ** 2))
        avg_m2_complex_error = np.nanmean(complex_error[:, 0])
        avg_complex_error_all_const = np.nanmean(complex_error_all_const)

        # wrap around 360 for phase
        for i, [ob, md] in enumerate(zip(ha_obs_results, ha_mod_results)):
            if not valid_idx[i]:
                continue
            for j, cons in enumerate(common_tidal_constituents):
                ob_idx = ob.name == cons
                md_idx = md.name == cons
                if abs(ob.g[ob_idx] - md.g[md_idx]) > 180:
                    if md.g[md_idx] < ob.g[ob_idx]:
                        md.g[md_idx] = md.g[md_idx] + 360
                    else:
                        md.g[md_idx] = md.g[md_idx] - 360

        with open(cache_file, 'wb') as f:
            pickle.dump(
                [ha_obs_results, ha_mod_results, stations, valid_idx, complex_error, complex_error_all_const],
                f, protocol=pickle.HIGHEST_PROTOCOL
            )
    else:
        with open(cache_file, 'rb') as f:
            ha_obs_results, ha_mod_results, stations, valid_idx, complex_error, complex_error_all_const = pickle.load(f)

    # plot 1: horizontally spread out stations
    fig = plt.figure(figsize=(30, 15))  # Adjust figsize as needed
    gs = gridspec.GridSpec(8, 1)
    n_subplot = 0
    for constit in ['M2', 'K1']:
        if constit == 'Max':
            ax = fig.add_subplot(gs[n_subplot:n_subplot + 2])
            ax.set_title('Max amplitude comparison')

            y_obs = [result.A[0] for result in ha_obs_results[valid_idx]]
            ax.plot(stations[valid_idx], y_obs, 'ro', label='Obs')
            y_mod = [result.A[0] for result in ha_mod_results[valid_idx]]
            ax.plot(stations[valid_idx], y_mod, 'b+', label='Mod')
            ax.legend()
            ax.tick_params(axis='x', rotation=90)
            ax.set_ylabel('Max amplitude (m)')
            ax.grid(True)
            n_subplot += 2
        else:
            ax = fig.add_subplot(gs[n_subplot:n_subplot + 2])
            ax.set_title(f'{constit} amplitude comparison')
            y_obs = [result.A[result.name == constit] for result in ha_obs_results[valid_idx]]
            ax.plot(stations[valid_idx], y_obs, 'ro', label='Obs')
            y_mod = [result.A[result.name == constit] for result in ha_mod_results[valid_idx]]
            ax.plot(stations[valid_idx], y_mod, 'b+', label='Mod')
            ax.legend()
            ax.tick_params(axis='x', rotation=90)
            ax.set_ylabel(f'{constit} amplitude (m)')
            ax.grid(True)
            n_subplot += 2

            ax = fig.add_subplot(gs[n_subplot])
            ax.set_title(f'{constit} phase comparison')
            y_obs = [result.g[result.name == constit] for result in ha_obs_results[valid_idx]]
            ax.plot(stations[valid_idx], y_obs, 'ro', label='Obs')
            y_mod = [result.g[result.name == constit] for result in ha_mod_results[valid_idx]]
            ax.plot(stations[valid_idx], y_mod, 'b+', label='Mod')
            ax.legend()
            ax.tick_params(axis='x', rotation=90)
            ax.set_ylabel(f'{constit} phase (deg)')
            ax.grid(True)
            n_subplot += 1

    # Plot complex error for 8 common tidal constituents at each station
    ax = fig.add_subplot(gs[n_subplot:])  # occupy the rest of the space
    ax.set_title(
        'Station-averaged complex error (8 consitituents):'
        f'{np.nanmean(complex_error_all_const):.3f} m'
    )
    # cumulative bar plot
    color = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange']
    bottom = np.zeros(len(stations[valid_idx]))
    for i, cons in enumerate(common_tidal_constituents):
        ax.bar(
            stations[valid_idx], complex_error[valid_idx, i], bottom=bottom, color=color[i],
            label=f'{cons}: {np.nanmean(complex_error[valid_idx, i]):.3f} m'
        )
        bottom += complex_error[valid_idx, i]
    
    ax.set_xlabel('Stations')
    ax.set_ylabel('Complex error (m)')
    ax.legend()
    ax.tick_params(axis='x', rotation=90)
    ax.grid(True)

    # Adjust layout
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(f'{event["cdir"]}/ha_{event["plot_start_day_str"]}_{event["plot_end_day_str"]}.svg')
    plt.show()

    # plot 2: scatter plot of model vs obs
    # 4 subplots: M2, K1, Amplitude, Phase
    # set global font size
    plt.rcParams.update({'font.size': 20})

    fig = plt.figure(figsize=(15, 15))  # Adjust figsize as needed
    n = 0  # number of subplots
    plot_consitituents = ['M2', 'K1']
    colors = [ (1, 0.7, 0.0), 'blue', 'red', 'black']
    ha_obs_valid_results = ha_obs_results[valid_idx]
    ha_mod_valid_results = ha_mod_results[valid_idx]
    valid_color_idx = color_idx[valid_idx]
    alpha = 0.5
    # Amplitude
    for constit in plot_consitituents:
        ax = fig.add_subplot(2, 2, n + 1)
        ax.set_title(f'{constit} amplitude (m) comparison')
        ax.set_xlabel('Obs')
        ax.set_ylabel('Mod')

        for i, group in enumerate(color_groups):
            mask = valid_color_idx == i
            y_obs = [result.A[result.name == constit] for result in ha_obs_valid_results[mask]]
            y_mod = [result.A[result.name == constit] for result in ha_mod_valid_results[mask]]
            ax.plot(y_obs, y_mod, 'o', alpha=alpha, color=colors[i], label=group, markersize=7)
            
        y_obs = [result.A[result.name == constit] for result in ha_obs_results[valid_idx]]
        ax.plot(y_obs, y_obs, 'k-', alpha=0.3)  # 45 degree line
        ax.grid(True)
        ax.legend()
        ax.set_aspect('equal', adjustable='box')
        n += 1

    # Phase
    for constit in plot_consitituents:
        ax = fig.add_subplot(2, 2, n + 1)
        ax.set_title(f'{constit} phase (degrees) comparison')
        ax.set_xlabel('Obs')
        ax.set_ylabel('Mod')

        for i, group in enumerate(color_groups):
            mask = valid_color_idx == i
            y_obs = [result.g[result.name == constit] for result in ha_obs_valid_results[mask]]
            y_mod = [result.g[result.name == constit] for result in ha_mod_valid_results[mask]]
            ax.plot(y_obs, y_mod, 'o', alpha=alpha, color=colors[i], label=group, markersize=7)

        y_obs = [result.g[result.name == constit] for result in ha_obs_results[valid_idx]]
        ax.plot(y_obs, y_obs, 'k-', alpha=0.3)  # 45 degree line
        ax.legend()
        ax.grid(True)
        ax.set_aspect('equal', adjustable='box')
        n += 1

    # Adjust layout
    plt.savefig(
        f'{event["cdir"]}/ha_scatter_'
        f'{event["plot_start_day_str"]}_{event["plot_end_day_str"]}.png'
    )
    plt.show()
    
    # plot 2.5: lower values in M2 amplitude
    plt.rcParams.update({'font.size': 30})
    constit = 'M2'
    fig = plt.figure(figsize=(8, 8))  # Adjust figsize as needed
    ax = fig.add_subplot(1, 1, 1)
    # ax.set_title(f'{constit} amplitude (m) comparison')
    ax.set_xlabel('Obs')
    ax.set_ylabel('Mod')

    ha_obs_valid_results = ha_obs_results[valid_idx]
    ha_mod_valid_results = ha_mod_results[valid_idx]
    valid_color_idx = color_idx[valid_idx]
    for i, group in enumerate(color_groups):
        mask = valid_color_idx == i
        y_obs = [result.A[result.name == constit] for result in ha_obs_valid_results[mask]]
        y_mod = [result.A[result.name == constit] for result in ha_mod_valid_results[mask]]
        ax.plot(y_obs, y_mod, 'o', alpha=alpha, color=colors[i], label=group, markersize=15)
            
    y_obs = [result.A[result.name == constit] for result in ha_obs_results[valid_idx]]
    ax.plot(y_obs, y_obs, 'k-', alpha=0.3)  # 45 degree line

    # limit x and y axis to 0.5
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 0.5)
    ax.grid(True)
    # Adjust layout
    plt.tight_layout()
    plt.savefig(
        f'{event["cdir"]}/ha_scatter_'
        f'{event["plot_start_day_str"]}_{event["plot_end_day_str"]}_zoomed_in.png'
    )
    plt.show()
    

if __name__ == '__main__':
    project_json = './Plot/ha.json'
    event_name = 'EX62'

    with open(project_json, 'r', encoding='utf-8') as f:
        main_dict = json.load(f)
    event = main_dict[event_name]
    
    # group stations by region
    station_groups = get_stations_from_groups(
        grouping_polygon_shpfile='/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/station_group_polygons.shp',
        station_bp_file=event['station_bp_file']
    )

    plot_ha(event, color_groups=station_groups, cache=True)
    print("done!")
