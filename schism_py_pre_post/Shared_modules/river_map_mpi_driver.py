import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from mpi4py import MPI
import gc
from schism_py_pre_post.Shared_modules.river_map_tif_preproc import find_thalweg_tile
from schism_py_pre_post.Shared_modules.make_river_map import make_river_map
from schism_py_pre_post.Grid.SMS import SMS_MAP
import time


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

'''
def plot_schism2D_parallel(
    rundir='./',
    model_start_time = datetime.strptime('2014-12-01 00:00:00', "%Y-%m-%d %H:%M:%S"),
    var_str=None, stacks=[1, 2, 3],
    caxis = [10, 35], output_dir='./'
):
    """
    function for plotting 2D variables of the schism outputs;
    mpi can be used when plotting many time stamps 
    """

    # devide tasks
    stacks_proc = stacks[my_mpi_idx(len(stacks), size, rank)]

    # stacks_proc = np.flip(stacks_proc)
    print(f'process {rank} handles {stacks_proc}\n')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    grd_name = f'{rundir}/hgrid.gr3'
    alt_grd_name = os.path.splitext(grd_name)[0] + '.pkl'
    if os.path.exists(alt_grd_name):
        gd = schism_grid(alt_grd_name)
    else:
        gd = schism_grid(grd_name)
        gd.save(alt_grd_name)

    for stack in stacks_proc:
        fname = f'{rundir}/outputs/{var_str}_{stack}.nc'
        print(f'Core {rank} reading {fname}\n')
        my_nc = netCDF4.Dataset(fname)
        var = my_nc.variables[var_str]
        this_time_dts = my_nc.variables['time']

        for j, this_time_dt in enumerate(this_time_dts):
            this_time = model_start_time + timedelta(seconds=this_time_dt.data.item())
            title_str = f"Surface {var_str} {this_time.strftime('%Y-%m-%d %H:%M:%S')}"
            savefilename = f'{output_dir}/{title_str.replace(" ", "_")}.png'

            # if os.path.exists(savefilename):
            #     print(f'Core {rank} skipping: {savefilename}\n')
            #     continue  # skip existing figures

            print(f'Core {rank} plotting time: {title_str} from {fname}\n\n')
            gd.dp = var[j, :, -1].filled(fill_value=-99999)
            gd.plot_grid(fmt=1, clim=caxis, levels=101, ticks=np.arange(caxis[0], caxis[1], 2), cmap='jet')

            plt.title(title_str)
            plt.gca().set_aspect('equal', 'box')
            plt.savefig(savefilename, dpi=400)
            # plt.show()
            plt.clf()
            plt.close()

        my_nc.close()
        del var, my_nc, model_start_time, this_time_dts
        gc.collect()

    print(f'Core {rank} finishing ...\n\n\n')
    comm.Barrier()
'''

if __name__ == "__main__":
    if rank == 0: print(f'A total of {size} cores used.')
    comm.barrier()

    dems_json_file = 'dems.json'
    thalweg_shp_fname='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/v13_thalweg_utm17n.shp'
    output_dir = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/SMS_MAP/'

    thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = \
        find_thalweg_tile(
            dems_json_file=dems_json_file,
            thalweg_shp_fname=thalweg_shp_fname,
            iNoPrint=bool(rank)  # only rank 0 prints to screen
        )
    if rank == 0: print(f'Thalwegs are divided into {len(tile_groups2thalwegs)} groups.')
    comm.barrier()

    # each core handles some groups
    my_idx = my_mpi_idx(N=len(tile_groups_files), size=size, rank=rank)
    my_tile_groups = tile_groups_files[my_idx]
    my_tile_groups_thalwegs = tile_groups2thalwegs[my_idx]
    print(f'Rank {rank} handles {np.squeeze(np.argwhere(my_idx))}')

    time_all_groups_start = time.time()

    for i, (my_tile_group, my_tile_group_thalwegs) in enumerate(zip(my_tile_groups, my_tile_groups_thalwegs)):
        time_this_group_start = time.time()
        try:
            make_river_map(
                tif_fnames = my_tile_group,
                thalweg_shp_fname = thalweg_shp_fname,
                thalweg_smooth_shp_fname = None,  # '/GA_riverstreams_cleaned_corrected_utm17N.shp'
                selected_thalweg = my_tile_group_thalwegs,
                output_dir = output_dir,
                output_prefix = f'{rank}_{i}_',
                mpi_print_prefix = f'[Rank {rank}, Group {i+1} of {len(my_tile_groups)}] ',
            )
        except:
            print(f'Rank {rank}: Group {i+1} of {len(my_tile_groups)} FAILED')
        
        print(f'Rank {rank}: Group {i+1} run time: {time.time()-time_this_group_start} seconds.')
    
    print(f'Rank {rank}: total run time: {time.time()-time_all_groups_start} seconds.')

    comm.Barrier()

    # write
    if rank == 0:
        map_files = glob.glob(f'{output_dir}/*_river.map')
        map_objects = [SMS_MAP(filename=map_file) for map_file in map_files]

        total_map = map_objects[0]
        for map_object in map_objects[1:]:
            total_map += map_object
        total_map.writer(filename=f'{output_dir}/total_arcs.map')

        print(f'Total run time: {time.time()-time_all_groups_start} seconds.')
