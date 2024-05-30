import os
import shutil

def link_OUT():
    '''
    create symlinks for transferring with Globus
    '''
    os.makedirs('OUT', exist_ok=True)
    os.makedirs('OUT/UV', exist_ok=True)
    os.makedirs('OUT/ST', exist_ok=True)
    os.makedirs('OUT/OZ', exist_ok=True)
    os.makedirs('OUT/HOT', exist_ok=True)

    os.chdir('OUT/UV')
    os.system('ln -s ../../outputs/horizontal*.nc .')

    os.chdir('../ST')
    os.system('ln -s ../../outputs/salinity_*.nc .')
    os.system('ln -s ../../outputs/temperature_*.nc .')

    os.chdir('../OZ')
    os.system('ln -s ../../outputs/out2d_*.nc .')
    os.system('ln -s ../../outputs/z*.nc .')

    os.chdir('../HOT')
    os.system('ln -s ../../outputs/hotstart_it*.nc .')

def link_sflux(wdir='./'):
    from glob import glob
    import os

    product_dict = {
        "hrrr": {
            "folder_patterns": [
                "[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]",  # date pattern, daily product
                "HRRR",  # combined product
            ],
        },
        "gfs": {
            "folder_patterns": [
                "[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]",  # date pattern, daily product
                "GFS",  # combined product
            ],
        }
    }

    products = ["gfs", "hrrr"]  # order matters
    vars = ['air', 'prc', 'rad']

    os.chdir(wdir)
    os.makedirs("sflux", exist_ok=True)
    os.chdir("sflux")

    for k, product in enumerate(products):
        for folder_pattern in product_dict[product]["folder_patterns"]:  # try different folder patterns
            files = sorted(glob(f"../{folder_pattern}/{product}*.nc"))
            if files:
                break
        if not files:
            raise FileNotFoundError(f"No files found for {product}")

        for i, file in enumerate(files):
            for var in vars:
                os.symlink(file, f"sflux_{var}_{k+1}.{i+1:04d}.nc")

    with open("sflux_inputs.txt", "w") as f:
        f.write("&sflux_inputs\n/\n")

def combine_hot(wdir):
    os.chdir(wdir)
    script = '/sciclone/home/feiye/bin/combine_hotstart7.viz'
    steps = [18432, 19584, 20736,]

    for step in steps:
        print(f'Processing hotstart at step {step}')
        os.system(f'{script} -i {step}') 


if __name__ == "__main__":
    combine_hot(wdir='/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15a/outputs/')