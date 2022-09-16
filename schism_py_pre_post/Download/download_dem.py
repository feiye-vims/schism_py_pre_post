import subprocess
import pathlib
import os


download_dir = r'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CuDEM/Lonlat/'

with open(f'{download_dir}/urllist8483.txt') as f:
    for line in f:
        line = line.rstrip()
        if pathlib.Path(line).suffix == ".tif":
            fname = os.path.basename(line)
            if ~os.path.exists('f{download_dir}/{fname}'):
                print(f'downloading {line}')
                process =subprocess.Popen(f'curl -O {line}', shell=True, cwd=download_dir)
                process.wait()
            else:
                print(f'{line} exists, skipping ...')
        pass