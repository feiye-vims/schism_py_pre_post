from schism_py_pre_post.Plot.plot_elev import plot_elev, get_hindcast_elev, get_obs_elev
from schism_py_pre_post.Grid.Bpfile import Bpfile
import pandas as pd
import json
import numpy as np
import os
from mpl_toolkits.basemap import Basemap  # mamba install basemap
import matplotlib.pyplot as plt


def pacific_stations(station_bp_file=None):
    # --------------define stations----------------------
    # stations, ICOGS v2 and v3, coastal act
    if station_bp_file is None:
        station_bp_file = '/sciclone/schism10/feiye/Coastal_Act/RUN17a/station.in'

    noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id

    stations_groups = {
        'Hawaii': noaa_stations_all[76:88],
        'West_Coast_1': noaa_stations_all[88:108],
        'West_Coast_2': noaa_stations_all[108:128],
        'West_Coast_3': noaa_stations_all[128:147],
    }
    default_datums = {
        'Hawaii': 'MSL',
        'West_Coast_1': 'NAVD',
        'West_Coast_2': 'NAVD',
        'West_Coast_3': 'NAVD',
    }
    return [stations_groups, default_datums]

def ecgc_stations(station_bp_file=None):
    # --------------define stations----------------------
    # stations, ICOGS v2 and v3, coastal act
    if station_bp_file is None:
        station_bp_file = '/sciclone/schism10/feiye/Coastal_Act/Station_in/15a_station.in'

    noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id

    stations_groups = {
        'Florida': noaa_stations_all[:10],
        'Atlantic': noaa_stations_all[10:29],
        'GoME': noaa_stations_all[29:39], 'GoMX_west': noaa_stations_all[40:60],
        'GoMX_east': noaa_stations_all[60:80], 'Atlantic_inland1': noaa_stations_all[80:100],
        'Atlantic_inland2': noaa_stations_all[100:120], 'GoMX_inland': noaa_stations_all[120:150],
        'Puerto_Rico': noaa_stations_all[150:163]
    }
    default_datums = {
        'Florida': 'NAVD',
        'Atlantic': 'NAVD',
        'GoME': 'NAVD', 'GoMX_west': 'NAVD',
        'GoMX_east': 'NAVD', 'Atlantic_inland1': 'NAVD',
        'Atlantic_inland2': 'NAVD', 'GoMX_inland': 'NAVD',
        'Puerto_Rico': 'MSL'
    }
    return [stations_groups, default_datums]


def subset_stations_in_box(box, station_bp_file, group_name="Landfall_region", default_datum='NAVD', max_subplots = 20):
    # select a subset of stations within the box and redefine noaa_stations_groups
    stations_all = Bpfile(station_bp_file, cols=5).st_id
    stations_lonlat = Bpfile(station_bp_file, cols=5).nodes[:, :2]
    idx = (
        (box['W']-1 < stations_lonlat[:, 0]) *
        (box['E']+1 > stations_lonlat[:, 0]) *
        (box['S']-1 < stations_lonlat[:, 1]) *
        (box['N']+1 > stations_lonlat[:, 1])
    )

    stations_inbox = np.array(stations_all)[idx]
    stations_groups = {}
    default_datums = {}
    for i_group, i in enumerate(range(0, len(stations_inbox), max_subplots)):
        these_stations = stations_inbox[i:i + max_subplots]
        this_group_name = f"{group_name}_{i_group}"
        stations_groups[this_group_name] = these_stations
        default_datums[this_group_name] = default_datum
        print(f'{len(these_stations)} stations in {this_group_name}')

    return [stations_groups, default_datums]


def stats_scatter(stats=None, stats_file=None, var_str=None,
                  region="Full_domain", box={"W": -100, "E": -60, "S": 8, "N": 48},
                  plot_symbol_dict=None, filename=None):
    if stats is None:
        stats = pd.read_csv(stats_file)
    var_str0 = var_str.split('_')[0]

    # plot stats as scatter points
    ilabel = plot_symbol_dict[region]["ilabel"]
    grid_spacing = plot_symbol_dict[region]["grid_spacing"]
    symbol_size = plot_symbol_dict[region]["symbol_size"]
    var_colorbar_str = plot_symbol_dict[var_str]["var_colorbar_str"]
    cmap = plot_symbol_dict[var_str]["cmap"]
    vmin = plot_symbol_dict[var_str]["vmin"]
    vmax = plot_symbol_dict[var_str]["vmax"]

    plt.figure(figsize=(12, 10), dpi=300)
    m = Basemap(projection="merc", resolution="f",  # mamba install basemap-data-hires
                llcrnrlat=box['S']-1, llcrnrlon=box['W']-1,
                urcrnrlat=box['N']+1, urcrnrlon=box['E']+1)
    m.shadedrelief()
    m.drawcoastlines()
    parallels = np.arange(np.floor(box['S'])-1, np.ceil(box['N'])+1, grid_spacing)
    m.drawparallels(parallels, labels=[False, True, True, False])
    meridians = np.arange(np.ceil(box['W'])-1, np.ceil(box['E'])+1, grid_spacing)
    m.drawmeridians(meridians, labels=[False, True, True, False])
    mx, my = m(np.array(stats['station_lon']), np.array(stats['station_lat']))    # transform coordinates
    plt.scatter(mx, my, symbol_size, stats[var_str0].values, cmap=cmap, vmin=vmin, vmax=vmax)
    if ilabel == 1:
        for i, label in enumerate(stats['station_id']):
            plt.annotate(f'{label}', (mx[i], my[i]))
    clb = plt.colorbar()
    clb.ax.set_title(f'{var_colorbar_str} (m)')
    plt.savefig(f'{filename}_{var_str}.png', dpi=400)
    plt.close('all')


if __name__ == "__main__":
    # ---------------------------------------------------------------------------------
    #    inputs
    # ---------------------------------------------------------------------------------
    hurricanes = ['Ida_v5']
    main_dict = 'stofs3d.json'  # 'coastal_act_stats_period_3D_1st_round.json'

    region = "Full_domain"  # "Landfall_region", "Full_domain", "Manual"
    var_str = 'MAE'
    nday_moving_average = 0

    with open(main_dict) as d:
        hurricane_dict = json.load(d)
    station_bp_file = hurricane_dict['All']['station_bp_file']

    other_dicts_files = []  # ['coastal_act_stats_period_3D_others.json']
    other_line_styles = ['g']
    other_shifts = [0]

    datum = ''  # '' (use predefined ones in ecgc_stations), 'NAVD', 'MSL'
    outfilename_suffix = 'Mostly_NAVD'  # 'Mostly_NAVD': note some stations don't have NAVD datum and their elevation time series will be demeaned and shown as "MSL"

    with open('coastal_act_stats_plot_symbols.json') as d:
        plot_symbol_dict = json.load(d)
    # --end inputs-------------------------------------------------------------------------------

    for hurricane in hurricanes:
        dict_item = hurricane_dict[hurricane]
        plot_start_day_str = dict_item['plot_start_day_str']
        plot_end_day_str = dict_item['plot_end_day_str']
        model_start_day_str = dict_item['model_start_day_str']
        elev_out_file = dict_item['elev_out_file']
        cdir = dict_item['cdir']
        runid = os.path.basename(os.path.dirname(cdir))
        if region == "Full_domain":
            box = hurricane_dict['All']['box']
            stations_groups, default_datums = ecgc_stations(station_bp_file)
            # stations_groups, default_datums = pacific_stations(station_bp_file)
        elif region == "Manual":
            stations_groups = hurricane_dict[hurricane]['stations_group']
            default_datums = hurricane_dict[hurricane]['default_datums']
        else:
            box = hurricane_dict[hurricane]['box']
            stations_groups, default_datums = subset_stations_in_box(box, station_bp_file, group_name=region, default_datum=datum)

        # ---------------------------------------------------------------------------------
        # # Manually overwrite the parameters read from '*.json'
        # overwrite default_datums
        if datum is not None and datum != '':
            for key in default_datums:
                default_datums[key] = datum

        # plot_start_day_str = '2018-08-17 00:00:00'
        # plot_end_day_str = '2018-10-01 00:00:00'
        # model_start_day_str = '2017-08-04 00:00:00'
        # # model outputs
        # elev_out_file = '/sciclone/schsm10/feiye/Coastal_Act/RUN13b/PostP/staout_1'
        # cdir = '$cdir/srv/www/htdocs/yinglong/feiye/Coastal_Act/RUN13b/'

        other_dicts = []
        for other_dict_file in other_dicts_files:
            with open(other_dict_file) as d:
                other_dicts.append(json.load(d))
        other_runids = [os.path.basename(os.path.dirname(x[hurricane]['cdir'])) for x in other_dicts]
        # ---------------------------------------------------------------------------------

        stats = pd.DataFrame()
        for group_name in stations_groups:
            filename_base = f'{hurricane}_{group_name}_{default_datums[group_name]}'
            # if group_name != 'Puerto_Rico':
            #     continue
            print(f'\n\n\n processing {group_name} for {hurricane}')
            # SCHISM's staout_1
            mod = get_hindcast_elev(
                model_start_day_str=model_start_day_str,
                noaa_stations=None,
                station_in_file=station_bp_file,
                elev_out_file=elev_out_file
            )

            # get obs
            [obs, datums, st_info] = get_obs_elev(
                plot_start_day_str=plot_start_day_str,
                plot_end_day_str=plot_end_day_str,
                noaa_stations=stations_groups[group_name],
                default_datum=default_datums[group_name]
            )

            # plot time series
            stat, fig_ax = plot_elev(obs, mod, plot_start_day_str, plot_end_day_str,
                                     stations_groups[group_name],
                                     datums, st_info, 'ts_' + filename_base, iplot=False, subplots_shape=(None, None),
                                     nday_moving_average=nday_moving_average, label_strs=['obs', runid])
            stats = stats.append(stat)
            if len(other_dicts) > 0:
                for other_runid, other_dict, other_line_style, other_shift in zip(other_runids, other_dicts, other_line_styles, other_shifts):
                    other_mod = get_hindcast_elev(
                        model_start_day_str=other_dict[hurricane]['model_start_day_str'],
                        noaa_stations=None,
                        station_in_file=station_bp_file,
                        elev_out_file=other_dict[hurricane]['elev_out_file']
                    )
                    plot_elev(obs, other_mod, plot_start_day_str, plot_end_day_str,
                              stations_groups[group_name],
                              datums, st_info, None, iplot=False, nday_moving_average=nday_moving_average,
                              fig_ax=fig_ax, line_styles=[None, other_line_style], shift=other_shift, label_strs=['obs', other_runid])
                    fig_ax[0].savefig(f'compare_ts_{filename_base}.png')

        # ---------------------------------------------------------------------------------
        filename_base = f'{hurricane}_{region}_{outfilename_suffix}'
        stats_scatter(stats=stats, var_str=var_str, region=region, plot_symbol_dict=plot_symbol_dict, filename=filename_base)

        rearranged_cols = ['station_lon', 'station_lat', 'station_id', 'RMSE', 'MAE',
                           'Bias', 'CC', 'Max_Obs', 'Max_Mod', 'station_name']
        stats = stats[rearranged_cols]
        stats.loc['mean'] = stats.iloc[:, :-1].mean()
        stats.at['mean', 'station_id'] = 'all'
        stats.at['mean', 'station_name'] = 'all'

        stats.to_csv(f'stats_{filename_base}.txt', encoding='utf-8', index=False, na_rep='-9999')

        # ---------------------------------------------------------------------------------

        # upload to ccrm drive:
        print(f"scp stats_{filename_base}.txt *png {cdir}/")
        os.system(f"scp stats_{filename_base}.txt *png {cdir}/")
        os.system(f"rm stats_{filename_base}.txt *png")
