from pylib import schism_grid
from hotstart_proc import Hotstart
import os
import numpy as np


if __name__ == "__main__":
    # input section
    wdir = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/'
    hycom_hot_file = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/hotstart.nc.0'
    my_hot_file = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/hotstart.nc'
    # end input section
    
    if os.path.exists(my_hot_file):
        os.system(f"rm {my_hot_file}")
    os.system(f"cp {hycom_hot_file} {my_hot_file}")
    
    my_hot = Hotstart(
        grid_info=wdir,
        hot_file=my_hot_file
    )

    for i, var in enumerate(['tem', 'sal']):
        hg = schism_grid(f'{wdir}/ecgc_coastal_{var}.gr3')
        idx = hg.dp > -9998
        for k in range(my_hot.grid.vgrid.nvrt):
            my_hot.tr_nd.val[idx, k, i] = hg.dp[idx]
    
    # set salinity to 0 for higher grounds
    rat=np.maximum(np.minimum(1.0, (my_hot.grid.hgrid.dp+3.0)/3.0), 0.0)
    hg.dp = rat
    hg.save('hgrid.gr3')
    my_hot.tr_nd.val[:, :, 1] *= np.transpose(np.tile(rat, (my_hot.grid.vgrid.nvrt, 1)))

    my_hot.writer(fname=my_hot_file)
