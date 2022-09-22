import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from mpi4py import MPI
import gc
from schism_py_pre_post.Rivermap.river_map_tif_preproc import find_thalweg_tile
from schism_py_pre_post.Rivermap.make_river_map import make_river_map
from schism_py_pre_post.Grid.SMS import SMS_MAP
import time
import pickle


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# print(f'process {rank} of {size}\n')

def my_mpi_idx(N, size, rank):
    my_idx = np.zeros((N, ), dtype=bool)
    n_per_rank, _ = divmod(N, size)
    n_per_rank = n_per_rank + 1
    # return slice())
    my_idx[rank*n_per_rank:min((rank+1)*n_per_rank, N)] = True
    return my_idx

if __name__ == "__main__":
    '''
    Driver for parallel execution of make_river_map.py.

    Thalwegs are grouped based on the DEM tiles associated with each thalweg.
    For each thalweg, its associated DEM tiles are those needed for determining
    the elevations on all thalweg points, as well as
    the elevations within a buffer zone of the thalweg (within which left and right banks will be sought)

    One core can be responsible for one or more thalweg groups,
    which are fed to make_river_map.py one at a time
    '''

    # ------------------------- input section --------------------------- 
    # files and dirs
    dems_json_file = 'dems.json'  # files for all DEM tiles

    thalweg_shp_fname='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/Shp/NWM_cleaned_ll_redist7m.shp'

    output_dir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/Outputs/' + \
        f'{os.path.basename(thalweg_shp_fname).split(".")[0]}_{size}cores/'

    # parameters 
    thalweg_buffer = 1000  # meters. This is the search range on either side of the thalweg.
                            # Because banks will be searched within this range,
                            # its value is needed now to associate DEM tiles with each thalweg
    cache_folder = '/sciclone/schism10/feiye/Cache/'
    i_DEM_cache = True  # Whether or not read DEM info from cache.
                        # The cache file saves box info of all DEM tiles.
                        # Reading from original tiles can be slow, so the default option is True
    i_thalweg_cache = True  # Whether or not read thalweg info from cache.
                             # The cache file saves coordinates, index, curvature, and direction at all thalweg points
                             # This is usually fast even without cache
    i_grouping_cache = True  # Whether or not read grouping info from cache,
                              # which is useful when the same DEMs and thalweg_shp_fname are used.
                              # A cache file named "dems_json_file + thalweg_shp_fname_grouping.cache"
                              # will be saved regardless of the option value.
                             # This is usually fast even without cache
    # ------------------------- end input section --------------------------- 

    if rank == 0: print(f'A total of {size} core(s) used.')
    comm.barrier()

    if rank == 0:
        if os.path.exists(output_dir):
            os.system(f"rm -r {output_dir}")
        else:
            os.mkdir(output_dir)

    # group thalwegs (with the option of cache)
    cache_name = cache_folder + \
        os.path.basename(dems_json_file) + '_' + \
        os.path.basename(thalweg_shp_fname) + '_grouping.cache'
    # Remove the existing cache if i_grouping_cache is False,
    # this assumes the cache file needs to be updated
    if (not i_grouping_cache) and os.path.exists(cache_name): os.remove(cache_name)
    # core 1 calculates the grouping, then saves a cache if the cache file does not exist
    if rank == 0:
        if not os.path.exists(cache_name):
            thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = \
                find_thalweg_tile(
                    dems_json_file=dems_json_file,
                    thalweg_shp_fname=thalweg_shp_fname,
                    iNoPrint=bool(rank), # only rank 0 prints to screen
                    cache_folder=cache_folder,
                    i_DEM_cache=i_DEM_cache, i_thalweg_cache=i_thalweg_cache
                )
            with open(cache_name, 'wb') as file:
                pickle.dump([thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs], file)
    comm.barrier()
    # all cores read from cache
    with open(cache_name, 'rb') as file:
        print(f'Reading grouping info from cache ...')
        thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = pickle.load(file)
    if rank == 0:
        print(f'Thalwegs are divided into {len(tile_groups2thalwegs)} groups.')
        for i, tile_group2thalwegs in enumerate(tile_groups2thalwegs):
            print(f'[ Group {i+1} ]-----------------------------------------------------------------------\n' + \
                  f'Group {i+1} includes the following thalwegs (idx starts from 0): {tile_group2thalwegs}\n' + \
                  f'Group {i+1} needs the following DEMs: {tile_groups_files[i]}\n')
    comm.barrier()

    # each core handles some groups
    my_idx = my_mpi_idx(N=len(tile_groups_files), size=size, rank=rank)
    # my_idx.fill(False); my_idx[[50]] = True
    my_tile_groups = tile_groups_files[my_idx]
    my_tile_groups_thalwegs = tile_groups2thalwegs[my_idx]
    print(f'Rank {rank} handles Group {np.squeeze(np.argwhere(my_idx))}\n')
    comm.Barrier()

    time_all_groups_start = time.time()

    for i, (my_tile_group, my_tile_group_thalwegs) in enumerate(zip(my_tile_groups, my_tile_groups_thalwegs)):
        time_this_group_start = time.time()
        # try:
        make_river_map(
            tif_fnames = my_tile_group,
            thalweg_shp_fname = thalweg_shp_fname,
            thalweg_smooth_shp_fname = None,  # '/GA_riverstreams_cleaned_corrected_utm17N.shp'
            selected_thalweg = my_tile_group_thalwegs,
            output_dir = output_dir,
            output_prefix = f'{rank}_{i}_',
            mpi_print_prefix = f'[Rank {rank}, Group {i+1} of {len(my_tile_groups)}] ',
        )
        # except:
        #     print(f'Rank {rank}: Group {i+1} of {len(my_tile_groups)} FAILED')

        print(f'Rank {rank}: Group {i+1} run time: {time.time()-time_this_group_start} seconds.')

    print(f'Rank {rank}: total run time: {time.time()-time_all_groups_start} seconds.')

    comm.Barrier()

    # write
    if rank == 0:
        map_files = glob.glob(f'{output_dir}/*_river.map')

        if len(map_files) > 0:
            map_objects = [SMS_MAP(filename=map_file) for map_file in map_files]

            total_map = map_objects[0]
            for map_object in map_objects[1:]:
                total_map += map_object
            total_map.writer(filename=f'{output_dir}/total_arcs.map')
        else:
            print('No map files found in final combination ...')

        print(f'Total run time: {time.time()-time_all_groups_start} seconds.')
