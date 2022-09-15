# %%    
import json
from statistics import median
from osgeo import gdal
from glob import glob
import copy
from schism_py_pre_post.Shared_modules.make_river_map import Tif2XYZ, get_all_points_from_shp
import numpy as np
import os
import matplotlib.pyplot as plt
import shapefile
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid

def get_tif_boxes(tif_files:list):
    tif_box = []
    for i, tif_file in enumerate(tif_files):
        # print(f'reading tifs: {i+1} of {len(tif_files)}, {tif_file}')
        tif = Tif2XYZ(tif_file)
        tif_box.append([min(tif.x), min(tif.y), max(tif.x), max(tif.y)])
    return tif_box

def reproject_tifs(tif_files:list, dstSRS='EPSG:26918'):
    for tif_file in tif_files:
        tif_outfile = os.path.splitext(tif_file)[0] + '.utm17N.tif'
        g = gdal.Warp(tif_outfile, tif_file, dstSRS=dstSRS)
        g = None

def pts_in_box(pts, box):
    in_box = (pts[:, 0] >  box[0]) * (pts[:, 0] <= box[2]) * \
             (pts[:, 1] >  box[1]) * (pts[:, 1] <= box[3])
    return in_box 

def find_parent_box(pts, boxes):
    parent = -np.ones((len(pts), 1), dtype=int)
    for j, box in enumerate(boxes):
        in_box = pts_in_box(pts[:,:2], box)
        parent[in_box] = j
    # plt.hist(parent, bins=len(np.unique(parent)))
    # np.savetxt('thalweg_parent.xyz', np.c_[pts[:,:2], parent])
    return parent

def tile2dem_file(dem_dict, dem_order, tile_code):
    DEM_id, tile_id = int(tile_code.real), int(tile_code.imag)
    if tile_id != -1:
        return dem_dict[dem_order[DEM_id]]['file_list'][tile_id]
    else:
        return None

def find_thalweg_tile(
    dems_json_file='dems.json',
    thalweg_shp_fname='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_riverstreams_cleaned_utm17N.shp',
):
    '''
    Assign thalwegs to DEM tiles
    '''
    # read DEMs
    with open(dems_json_file) as d:
        dem_dict = json.load(d)

    # get the box of each tile of each DEM
    dem_order = []
    for k, v in dem_dict.items():
        dem_order.append(k)
        dem_dict[k]['file_list'] = glob(dem_dict[k]['glob_pattern'])
        dem_dict[k]['boxes'] = get_tif_boxes(dem_dict[k]['file_list'])

    # read thalwegs
    thalweg_shp_fname = thalweg_shp_fname
    xyz, l2g, curv = get_all_points_from_shp(thalweg_shp_fname)

    # find DEM tiles for all thalwegs' points
    thalwegs2dems = [find_parent_box(xyz[:,:2], dem_dict[k]['boxes']) for k in dem_dict.keys()]

    # find DEM tiles for each thalweg
    thalwegs = []
    thalwegs_parents = []
    for i, idx in enumerate(l2g):
        line = xyz[idx,:]
        thalwegs.append(line)

        thalweg_parents = []  # one thalweg can have parent tiles from all DEM sources
        for i_dem, thalwegs2dem in enumerate(thalwegs2dems):
            thalweg2dem = np.unique(thalwegs2dem[idx]).tolist()
            thalweg_parents += [complex(i_dem, x) for x in thalweg2dem]  # real part is DEM id; complex part is tile id

        thalwegs_parents.append(thalweg_parents)

    # Group thalwegs: thalwegs from the same group have the same parent tiles
    groups = []
    group_id = 0
    thalweg2group = -np.ones((len(thalwegs)), dtype=int)
    for i, thalweg_parents in enumerate(thalwegs_parents):
        if thalweg_parents not in groups:
            groups.append(thalweg_parents)
            thalweg2group[i] = group_id
            group_id += 1
        else:
            for j, x in enumerate(groups):
                if x == thalweg_parents:
                    thalweg2group[i] = j

    # reduce groups: merge smaller groups into larger groups
    groups = np.array(groups, dtype=object)
    ngroup = len(groups)
    grp2large_grp = ngroup * np.ones((ngroup+1,), dtype=int)  # add a dummy mapping at the end
    for i1, group1 in enumerate(groups):
        for i2, group2 in enumerate(groups):
            if len(group1) < len(group2) and all(elem in group2 for elem in group1):
                # print(f'{group1} is contained in {group2}')
                grp2large_grp[i1] = i2
                break

    # But some large groups are still contained in larger groups,
    # get to the bottom of the family tree (e.g., parent's parent's parent ...)
    parents = np.squeeze(grp2large_grp[grp2large_grp])  # parent's parent
    idx = parents != len(groups)  # where parent's parent exists
    while any(idx):
        grp2large_grp[idx] = parents[idx]  # reset parent to parent's parent
        parents = parents[parents]  # advance family tree
        idx = parents != len(groups)  # see where parent's parent still exists

    idx = grp2large_grp==len(groups)  # where parent's parent is no-existent
    grp2large_grp[idx] = np.arange(len(groups)+1)[idx]  # parent group is self
    grp2large_grp = grp2large_grp[:-1]  # remove the dummy group at the end

    large_groups = groups[np.unique(grp2large_grp)]
    print(f'number of groups after reduction: {len(large_groups)}')
    group_lens = [len(x) for x in large_groups]
    print(f'group lengths: min {min(group_lens)}; max {max(group_lens)}; mean {np.mean(group_lens)}')

    thalweg2group = grp2large_grp[thalweg2group]
    map_grp = dict(zip(np.unique(grp2large_grp), np.arange(len(np.unique(grp2large_grp)))))
    thalweg2large_group = np.array([map_grp[x] for x in thalweg2group])

    large_group2thalwegs = [[] for _ in range(len(large_groups))]
    for i, x in enumerate(thalweg2large_group):
        large_group2thalwegs[x].append(i)
    
    large_groups_files = copy.deepcopy(large_groups)
    for i, group in enumerate(large_groups):
        for j, tile_code in enumerate(group):
            large_groups_files[i][j] = tile2dem_file(dem_dict=dem_dict, dem_order=dem_order, tile_code=tile_code)
            
    # histogram
    # plt.hist(thalweg2large_group, bins=len(np.unique(thalweg2large_group)))
    # plt.show()

    return thalweg2large_group, large_groups_files, np.array(large_group2thalwegs, dtype=object)

if __name__ == "__main__":
    # find_thalweg_tile()
    # %%
    # Reproject
    # tif_files = glob(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CRM/Lonlat/*.tif')
    # reproject_tifs(tif_files, 'EPSG:26917')

    # Merge small coned tiles into larger ones (similar to CuDEM's tile size)
    # cudem_files = glob(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CuDEM/*.tif')
    # cudem_boxes = get_tif_boxes(cudem_files)

    # coned_files = glob(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CoNED/Original/*.tif')
    # coned_boxes = get_tif_boxes(coned_files)

    # coned_centers = np.c_[
    #     (np.array(coned_boxes)[:, 0]+np.array(coned_boxes)[:, 2])/2, 
    #     (np.array(coned_boxes)[:, 1]+np.array(coned_boxes)[:, 3])/2,
    # ]

    # for i, cudem_box in enumerate(cudem_boxes):
    #     in_box = (coned_centers[:, 0] >  cudem_box[0]) * \
    #              (coned_centers[:, 0] <= cudem_box[2]) * \
    #              (coned_centers[:, 1] >  cudem_box[1]) * \
    #              (coned_centers[:, 1] <= cudem_box[3])
    #     in_box_files = (np.array(coned_files)[in_box]).tolist()
    #     g = gdal.Warp(f"/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CoNED/Combined/GA_CoNED_merged_{i}.tif",
    #                   in_box_files, format="GTiff", options=["COMPRESS=LZW", "TILED=YES"])
    #     g = None # Close file and flush to disk

    pass

