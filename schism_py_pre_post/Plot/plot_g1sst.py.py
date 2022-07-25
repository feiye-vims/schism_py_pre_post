from schism_py_pre_post.Download import download_g1sst
from schism_py_pre_post.Datasets import G1SST
import glob
import os
import matplotlib.pyplot as plt
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def my_mpi_idx(N, size, rank):
    if N > size:
        n_per_rank, _ = divmod(N, size)
        n_per_rank = n_per_rank + 1
        return slice(rank*n_per_rank, min((rank+1)*n_per_rank, N))
    else:
        if rank+1 > N:
            return slice

def parallel_plot(rank, size) -> None:
    # -----------------------------------------------------------------------------------
    # ------------plot G1SST----------------
    # -----------------------------------------------------------------------------------
    # plot all days
    files = glob.glob(f'{wdir}/*.nc')
    my_idx = my_mpi_idx(len(files), size, rank)
    my_files = files[my_idx]

    for file in my_files:
        if os.path.exists(file[:-2] + "png"):
            continue
        single_plot_g1sst(file)

def single_plot_g1sst(fname) -> None:
    # plot a single day
    data = G1SST(fname)
    data.plot_2d(
        i_show_plot=1, this_var='analysed_sst',
        t_idx=0, this_time_str=None,
        xlim=[-98, -60], ylim=[8.5, 46], clim=[10, 35],
        model_basetime='2014-12-01T00:00:00Z',  # used for marking model days in the title
    )
    plt.savefig(fname[:-2] + 'png', dpi=600)
    plt.clf()
    plt.close()
    del data

if __name__ == "__main__":
    wdir = '/sciclone/data10/feiye/G1SST/2015/'

    # download_g1sst(download_dir=wdir, begin_date='2014-12-01', end_date='2016-01-01')

    # parallel_plot(rank, size)

    single_plot_g1sst(f'{wdir}/20160101.nc')
    pass
