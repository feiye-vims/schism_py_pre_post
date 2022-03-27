import os
import pandas as pd
from schism_py_pre_post.Grid.SMS import SMS_MAP, SMS_ARC
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
        my_arc = SMS_ARC(points=np.array(arc['geometry']['coordinates']))
        arc_list.append(my_arc)
        # my_arc_hats = my_arc.make_hats(arc_hat_length=30/110844.0420526655)
        # arc_hat_list += my_arc_hats

        xyz = np.r_[xyz, my_arc.points]

    return SMS_MAP(arcs=arc_list), xyz


if __name__ == "__main__":
    # Specify levee ids
    wdir = '/sciclone/scr10/feiye/NLD/'
    levee_name = 'FL_levees'
    df = pd.read_csv(f'{wdir}/{levee_name}_info.csv')
    system_ids = df['System_ID'].to_numpy().astype(int).tolist()  # system_ids = [4405000525, 1605995181, 5905000001]

    # Download profiles
    download_nld(output_dir=wdir, output_fname=levee_name, levee_id_list=system_ids)

    # convert to sms map
    my_map, xyz = nld2map(nld_fname=f'{wdir}/{levee_name}/System.geojson')
    my_map.writer(filename=f'{wdir}/{levee_name}/{levee_name}.map')
    # my_map_hats.writer('test_hats.map')

    pass 
