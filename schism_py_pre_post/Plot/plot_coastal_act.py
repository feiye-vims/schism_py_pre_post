from schism_py_pre_post.Download.download_coops_elev import get_coops_elev
from schism_py_pre_post.Plot.plot_elev import plot_elev, get_hindcast_elev, datum_shift
from schism_py_pre_post.Grid.Bpfile import Bpfile
import pandas as pd
import json
import numpy as np
import os
from mpl_toolkits.basemap import Basemap  # mamba install basemap
import matplotlib.pyplot as plt
from glob import glob


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

sample_stations = {'Sample': ['8720030', '8423898', '8419317', '8410140']}
sample_datum = {'Sample': 'NAVD'}


def ecgc_stations_subset(station_bp_file=None):
    # --------------define stations----------------------
    # stations, ICOGS v2 and v3, coastal act
    if station_bp_file is None:
        station_bp_file = '/sciclone/schism10/feiye/Coastal_Act/Station_in/15a_station.in'

    noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id

    stations_groups = {
        'Florida': noaa_stations_all[:10], 'GoME': noaa_stations_all[29:39],
        # 'Atlantic': noaa_stations_all[10:29], 'Atlantic_inland1': noaa_stations_all[80:100], 'Atlantic_inland2': noaa_stations_all[100:120],
        'GoMX_inland': noaa_stations_all[120:150], 'GoMX_east': noaa_stations_all[60:80], 'GoMX_west': noaa_stations_all[41:60],
        'Puerto_Rico': noaa_stations_all[150:164] + noaa_stations_all[39:41]  # last 2 are bermuda
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

def ecgc_stations(station_bp_file=None):
    # --------------define stations----------------------
    # stations, ICOGS v2 and v3, coastal act
    if station_bp_file is None:
        station_bp_file = '/sciclone/schism10/feiye/Coastal_Act/Station_in/15a_station.in'

    noaa_stations_all = Bpfile(station_bp_file, cols=5).st_id

    stations_groups = {
        'Florida': noaa_stations_all[:10],
        'Atlantic': noaa_stations_all[10:29],
        'GoME': noaa_stations_all[29:39], 'GoMX_west': noaa_stations_all[41:60],
        'GoMX_east': noaa_stations_all[60:80], 'Atlantic_inland1': noaa_stations_all[80:100],
        'Atlantic_inland2': noaa_stations_all[100:120], 'GoMX_inland': noaa_stations_all[120:150],
        'Puerto_Rico': noaa_stations_all[150:164] + noaa_stations_all[39:41]  # last 2 are bermuda
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

def write_stat(stats, fname):
    rearranged_cols = ['station_name', 'station_id', 'station_lon', 'station_lat', 'RMSE', 'MAE',
                    'Bias', 'CC', 'ubRMSE', 'Max_Obs', 'Max_Mod']
    stats = stats[rearranged_cols]
    means = stats.iloc[:, 2:].mean(axis=0)
    new_row = pd.Series(['Mean'] * len(stats.columns), index=stats.columns)
    new_row[2:] = means.values
    stats = pd.concat([pd.DataFrame([new_row]), stats], ignore_index=True)

    # stats.to_csv(fname, encoding='utf-8', index=False, na_rep='-9999', float_format='%.4f')
    pd.options.display.float_format = '{:,.4f}'.format
    formatted_string = stats.to_string(index=False)
    with open(fname, 'w') as f:
        f.write(formatted_string)

    return stats.iloc[:1, :].to_string(index=False).split('\n')  # return the first row (var names) and second row (values)

if __name__ == "__main__":
    # ---------------------------------------------------------------------------------
    #    inputs
    # ---------------------------------------------------------------------------------
    events = ['STOFS3D-v7_reforecast3']

    datum = 'NAVD'  # '': (use predefined ones in ecgc_stations); 'NAVD'; 'MSL'
    shift = 0.0  # uniform shift for all stations
    # if model outputs are not in NAVD, shift them to NAVD; None if no shift
    datum_shift_file = '/sciclone/schism10/feiye/STOFS3D-v6/fcst/run/navd2xgeoid_shift.txt' # None
    nday_moving_average = 0  # 2-day moving average

    region = "Full_domain"  # "Sample", "Landfall_region", "Full_domain", "FromDict"

    # neglect stations, only for stats
    # neglect_stations =['8551910', '8548989', '8545240', '8539094', '8720226', '8652857', '8447435', '8575512', '8738043', '8729108']
    # neglect_stations = [
    #     '8720030',  # 
    #     '8767816',  # Lake Charles
    #     '8423898', '8419317', '8410140'  # GoME
    # ]
    neglect_stations = []  # Washington DC, Lake Charles

    var_str = 'MAE'  # for the scatter plot

    main_dict = '/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/stofs3d.json'  # 'coastal_act_stats_period_3D_1st_round.json'
    other_dicts_files = []  # ['/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/stofs3d_other.json']  # ['coastal_act_stats_period_3D_others.json']
    other_shifts = []  # uniform shift for all stations and for each of the other runs
    other_line_styles = ['g']

    outfilename_suffix = ''  # 'Mostly_NAVD': some stations don't have NAVD datum and their elevation time series will be demeaned and shown as "MSL"

    subplots_shape = [10, None]  # n per col
    # --end inputs-------------------------------------------------------------------------------

    with open('/sciclone/data10/feiye/schism_py_pre_post/schism_py_pre_post/Plot/coastal_act_stats_plot_symbols.json') as d:
        plot_symbol_dict = json.load(d)

    cache_folder = os.path.realpath(os.path.expanduser('~/schism10/Cache/'))

    with open(main_dict) as d:
        hurricane_dict = json.load(d)
    station_bp_file = hurricane_dict[events[0]]['station_bp_file']
    with open(station_bp_file) as f:
        f.readline()
        n_station = int(f.readline().split()[0])
    station_subset = range(n_station)
    other_subsets = [range(n_station)]

    for event in events:
        dict_item = hurricane_dict[event]
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
        elif region == "FromDict":
            stations_groups = hurricane_dict[event]['stations_group']
            default_datums = hurricane_dict[event]['default_datums']
            try:
                station_subset = hurricane_dict[event]['station_subset']
            except KeyError:
                print('No station_subset in the dict, assuming all stations in station.in are used')
            
        elif region == "Sample":
            stations_groups = sample_stations
            default_datums = sample_datum
        else:
            box = hurricane_dict[event]['box']
            stations_groups, default_datums = subset_stations_in_box(box, station_bp_file, group_name=region, default_datum=datum)

        # ---------------------------------------------------------------------------------
        # # Manually override the parameters read from '*.json', only for testing
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

        other_dicts = []; other_runs_stats = [pd.DataFrame()] * len(other_dicts_files)
        for other_dict_file in other_dicts_files:
            with open(other_dict_file) as d:
                other_dicts.append(json.load(d))
        other_runids = [os.path.basename(os.path.dirname(x[event]['cdir'])) for x in other_dicts]
        # ---------------------------------------------------------------------------------

        # get all obs and save cache, if the connection is unstable,
        # there will be exception, then rerun this script until all obs are retrieved
        nvalid = 0
        n_station = 0

        dry_run = False
        if dry_run:
            for iter_groups, [group_name, stations] in enumerate(stations_groups.items()):
                [obs, datums, st_info] = get_coops_elev(
                    plot_start_day_str, plot_end_day_str, stations,
                    default_datum=default_datums[group_name],
                    cache_folder=cache_folder,
                    retrieve_method='searvey'
                )
                n_station += len(stations)
                nvalid += sum([ob is not None for ob in obs])
            print(f'{nvalid} out of {n_station} stations have valid data')

        final_datums = []
        stats = pd.DataFrame()
        group_stats = ''
        other_group_stats = ['' for _ in other_runids]
        for iter_groups, [group_name, stations] in enumerate(stations_groups.items()):

            filename_base = f'{event}_{group_name}_{default_datums[group_name]}'

            # stations = [station for station in stations if station not in neglect_stations]

            # if group_name != 'Puerto_Rico':
            #     continue
            print(f'\n\n\n processing {group_name} for {event}')
            # SCHISM's staout_1
            mod = get_hindcast_elev(
                model_start_day_str=model_start_day_str,
                noaa_stations=None,
                station_in_file=station_bp_file,
                elev_out_file=elev_out_file,
                station_in_subset=station_subset,
            )
            
            mod += shift
            if datum_shift_file is not None:
                mod = datum_shift(mod, datum_shift_file=datum_shift_file)

            # get obs
            [obs, datums, st_info] = get_coops_elev(
                begin_time=plot_start_day_str,
                end_time=plot_end_day_str,
                noaa_stations=stations,
                default_datum=default_datums[group_name],
                cache_folder=cache_folder,
                retrieve_method='noaa_coops'
            )
            final_datums += datums

            # plot time series
            stat, fig_ax = plot_elev(obs, mod, plot_start_day_str, plot_end_day_str, stations,
                                     datums, st_info, 'ts_' + filename_base, iplot=False, subplots_shape=subplots_shape,
                                     nday_moving_average=nday_moving_average, label_strs=['obs', runid])
            stats = pd.concat([stats, stat], axis=0, ignore_index=True)
            mean_stats_string = write_stat(stat, f'stats_{runid}_{group_name}.txt')
            if iter_groups == 0:  # add header before the first group
                group_stats += ''.ljust(27) + mean_stats_string[0] + "\n"
            group_stats += f'{group_name.ljust(25)}: {mean_stats_string[1]}\n'

            if len(other_dicts) > 0:
                for i, [other_runid, other_dict, other_line_style, other_shift, other_subset] in enumerate(zip(other_runids, other_dicts, other_line_styles, other_shifts, other_subsets)):
                    other_mod = get_hindcast_elev(
                        model_start_day_str=other_dict[event]['model_start_day_str'],
                        noaa_stations=None,
                        station_in_file=station_bp_file,
                        elev_out_file=other_dict[event]['elev_out_file'],
                        station_in_subset=other_subset
                    )
                    other_mod += other_shift
                    other_mod = datum_shift(other_mod, datum_shift_file='/sciclone/schism10/feiye/STOFS3D-v6/fcst/run/navd2xgeoid_shift.txt')
                    other_stat, _ = plot_elev(obs, other_mod, plot_start_day_str, plot_end_day_str,
                                              stations_groups[group_name],
                                              datums, st_info, None, iplot=False, nday_moving_average=nday_moving_average, subplots_shape=subplots_shape,
                                              fig_ax=fig_ax, line_styles=[None, other_line_style], label_strs=['obs', other_runid])
                    other_runs_stats[i] = pd.concat([other_runs_stats[i], other_stat], axis=0, ignore_index=True)

                    mean_stats_string = write_stat(other_stat, f'stats_{other_runid}_{group_name}.txt')
                    if iter_groups == 0:  # add header before the first group
                        other_group_stats[i] +=  ''.ljust(25) + mean_stats_string[0] + "\n"
                    other_group_stats[i] += f'{group_name.ljust(25)}: {mean_stats_string[1]}\n'

                    fig_ax[0].savefig(f'compare_ts_{filename_base}.png')

        # ---------------------------------------------------------------------------------
        filename_base = f'{event}_{region}_{outfilename_suffix}'
        stats_scatter(stats=stats, var_str=var_str, region=region, plot_symbol_dict=plot_symbol_dict, filename=filename_base)

        # mask neglected stations
        mask = ~stats['station_id'].isin(neglect_stations)

        mean_stats_string = write_stat(stats[mask], f'stats_{runid}_{filename_base}.txt')
        group_stats += "Overall".ljust(25) + ":" + mean_stats_string[1]
        for i, [other_runid, other_run_stats] in enumerate(zip(other_runids, other_runs_stats)):
            mean_stats_string = write_stat(other_run_stats[mask], f'stats_{other_runid}_{filename_base}.txt')
            other_group_stats[i] += "Overall".ljust(25) + ":" + mean_stats_string[1]

        for this_run_id, stats_string in zip([runid] + other_runids, [group_stats] + other_group_stats):
            with open(f'mean_stats_{this_run_id}_{filename_base}.txt', 'w') as f:
                f.write(stats_string)
        
        overall_stats_datum_info_texts = [
            f"Stations with NAVD datum: {sum(np.array(final_datums)=='NAVD')}",
            f"Stations with MSL datum: {sum(np.array(final_datums)=='MSL')}",
            f"Stations without data: {sum(np.array(final_datums)==None)}",
        ]

        # fname = 'stats_datum_info.txt'
        # with open(fname, 'w') as f:
        #     os.system(f"head -n 1 stats_{runid}_{filename_base}.txt > {fname}")
        #     os.system(f"tail -n 1 stats_{runid}_{filename_base}.txt > {fname}")
        #     for group in groups:
        #         f.write('\n'.join(overall_stats_datum_info_texts))
        
        # ---------------------------------------------------------------------------------

        # pause for key press
        input("Press Enter to continue...")

        # upload to ccrm drive:
        os.system(f"scp *stats*txt *png {cdir}/")
        os.system(f"rm *stats*txt *png")
        
        pass
            