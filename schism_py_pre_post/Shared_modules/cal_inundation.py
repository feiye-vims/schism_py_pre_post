from pylib import schism_grid, inside_polygon  # from ZG's pylib: pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pylibs4schism==0.1.10
import matplotlib.pyplot as plt
import numpy as np
from time import time
from schism_py_pre_post.Geometry.inpoly import find_ele_node_in_shpfile
import copy
import pandas as pd


if __name__ == "__main__":
    runids = ['RUN21a', 'RUN21f', 'RUN21l', 'RUN21m']
    maxelev_fnames = [
        'elev_max_stack121.1-121.24.gr3', 'elev_max_stack121.1-121.24.gr3',
        'elev_max_stack6.1-6.24.gr3', 'elev_max_stack6.1-6.24.gr3',
    ]

    for runid, maxelev_fname in zip(runids, maxelev_fnames):
        wdir = f'/sciclone/schism10/feiye/ICOGS/{runid}/PostP/Maxelev/'
        shapefile_name = f'/sciclone/schism10/feiye/ICOGS/Shapefiles/Orleans_levee.shp'
        maxelev_fpathname = f'{wdir}/{maxelev_fname}'

        maxelev_gd = schism_grid(maxelev_fpathname)
        hgrid = schism_grid('/sciclone/schism10/feiye/ICOGS/RUN21a/hgrid.utm.pkl')

        hgrid.compute_area()
        ele_dp = hgrid.interp_node_to_elem()

        inun_depth = maxelev_gd.dp + hgrid.dp
        inun_gd = copy.deepcopy(hgrid)
        inun_gd.dp = inun_depth
        inun_ele = inun_gd.interp_node_to_elem()

        [inpoly_ele_idx_list, _] = find_ele_node_in_shpfile(shapefile_name=shapefile_name, grid=hgrid)
        inun_props = []
        for inpoly_ele_idx in inpoly_ele_idx_list:
            sum_area = sum(hgrid.area[inpoly_ele_idx])
            sig_inun = inun_ele > 0.3048
            inun_area = sum(hgrid.area[inpoly_ele_idx & sig_inun])
            inun_volume = sum(hgrid.area[inpoly_ele_idx & sig_inun]*inun_ele[inpoly_ele_idx & sig_inun])
            max_inun = max(inun_ele[inpoly_ele_idx & sig_inun])

            inun_prop = [sum_area, inun_area, inun_area/sum_area, inun_volume, inun_volume/inun_area, max_inun]
            inun_props.append(inun_prop)
        inun_props = np.array(inun_props)
        total_sum_area = sum(inun_props)

        df = pd.DataFrame(inun_props, columns=['Sum_area', 'Inun_area', 'Percent_inun_area', 'Inun_volume', 'Avg_inun_depth', 'Max_inun_depth'])
        df.loc['Total'] = df.sum()
        df['Percent_inun_area']['Total'] = df['Inun_area']['Total']/df['Sum_area']['Total']
        df['Avg_inun_depth']['Total'] = df['Inun_volume']['Total']/df['Inun_area']['Total']
        df['Max_inun_depth']['Total'] = max(df['Max_inun_depth'].iloc[:-1])

        df.to_csv(f"{maxelev_fpathname}.inun_prop.txt")
        pass
