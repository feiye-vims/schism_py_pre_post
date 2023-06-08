import xarray as xr
import os
from glob import glob


if __name__ == "__main__":
    sflux_dir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I13p2/sflux/'
    sflux_outdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I13p2/sflux+3/'

    sflux_set = 'air_2'
    sflux_files = sorted(glob(f'{sflux_dir}/sflux_{sflux_set}.*.nc'))

    for sflux_file in sflux_files:
        print(sflux_file)
        sflux = xr.open_dataset(sflux_file)
        sflux['stmp'] += 2
        sflux.to_netcdf(f'{sflux_outdir}/{os.path.basename(sflux_file)}')
        sflux.close()
    
    pass