import os
import pandas as pd
from schism_py_pre_post.Grid.SMS import SMS_MAP, SMS_ARC, Levee_SMS_MAP
import numpy as np


def download_nld(output_dir=None, output_fname=None, levee_id_list=[]):
    if output_dir is None:
        output_dir = './'
    if output_fname is None:
        output_fname = 'test'
    if levee_id_list == []:
        raise Exception('no levee IDs specified')
        
    full_name = f'{output_dir}/{output_fname}'
    if not os.path.exists(f'{full_name}.zip'):
        curl_cmd = ('curl -X POST "https://levees.sec.usace.army.mil:443/api-local/download/dataset/geojson.zip"'
                    ' -H "accept: application/json" -H "Content-Type: application/json"' 
                    f' -d "{levee_id_list}" --output {full_name}.zip')
        os.system(curl_cmd)
    if not os.path.isdir(full_name):
        os.system(f'unzip {full_name}.zip -d {full_name}')

def nld2map(nld_fname=None):
    data = pd.read_json(nld_fname)

    arc_list = []
    # arc_hat_list = []
    xyz = np.zeros((0, 3), dtype=float)
    for arc in data['features']:
        points = np.array(arc['geometry']['coordinates'])
        my_arc = SMS_ARC(points=points, src_prj='epsg:4326')
        arc_list.append(my_arc)

        xyz = np.r_[xyz, my_arc.points[:, :3]]

    return SMS_MAP(arcs=arc_list), xyz

if __name__ == "__main__":
    # Specify levee ids
    wdir = '/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23/Grids/'
    levee_name = 'TX_levees'
    df = pd.read_csv(f'{wdir}/{levee_name}_info.csv')
    system_ids = df['System_ID'].to_numpy().astype(int).tolist()  # system_ids = [4405000525, 1605995181, 5905000001]

    # Download profiles
    download_nld(output_dir=wdir, output_fname=levee_name, levee_id_list=system_ids)

    # convert to sms map
    my_map, xyz = nld2map(nld_fname=f'{wdir}/{levee_name}/System.geojson')
    # my_leveemap = Levee_SMS_MAP(arcs=my_map.arcs, epsg=4326)
    # scale = 110852.4248 * 285.0 / 300.0
    # centerline_map, offsetline_map = my_leveemap.make_levee_maps(subsample=[300/scale, None], offset_list=[-3/scale, 3/scale, -12/scale, 12/scale])

    my_map.writer(filename=f'{wdir}/{levee_name}/{levee_name}.map')
    # centerline_map.writer(filename=f'{wdir}/{levee_name}/{levee_name}_centerline_300m.map')
    # offsetline_map.writer(filename=f'{wdir}/{levee_name}/{levee_name}_offsetline_300m.map')

    # my_map_hats.writer('test_hats.map')

    pass 
