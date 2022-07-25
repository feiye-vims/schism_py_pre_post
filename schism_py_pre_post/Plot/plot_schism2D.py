from pylib import schism_grid
import netCDF4
from datetime import datetime, timedelta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from mpi4py import MPI
import gc


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print(f'process {rank} of {size}\n')

    
def my_mpi_idx(N, size, rank):
    if N > size:
        n_per_rank, _ = divmod(N, size)
        n_per_rank = n_per_rank + 1
        return slice(rank*n_per_rank, min((rank+1)*n_per_rank, N))
    else:
        if rank+1 > N:
            return slice

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

if __name__ == "__main__":
    rundir = '/sciclone/schism10/feiye/STOFS3D-v4/RUN23m/'
    model_start_time = datetime.strptime('2014-12-01', "%Y-%m-%d")
    var_str = 'temperature'
    stacks = np.arange(0, 397)
    caxis = [10, 35]
    output_dir = '/sciclone/schism10/feiye/STOFS3D-v4/Outputs/O23m/Plots/'

    plot_schism2D_parallel(rundir, model_start_time, var_str, stacks, caxis, output_dir)

