# %%
from pylib import schism_grid, sms2grd, grd2sms
import os
import numpy as np
from schism_py_pre_post.Grid.Bpfile import Bpfile
from schism_py_pre_post.Grid.Hgrid_ported import read_schism_hgrid_cached
from schism_py_pre_post.Grid.SMS import SMS_MAP, lonlat2cpp, cpp2lonlat
from schism_py_pre_post.Shared_modules.set_levee_profile import set_levee_profile
from schism_py_pre_post.Shared_modules.set_additional_dp import set_additional_dp_v11_91
import pathlib
import copy
import pickle

from sqlalchemy import over

# %%  
def find_large_small_dp():
    gd = sms2grd('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/maxelev.mask.2dm')
    sorted_dp = np.sort(gd.dp)
    sorted_idx = np.argsort(gd.dp)
    valid = np.squeeze(np.argwhere(~np.isnan(sorted_dp)))
    n = 30
    print(np.c_[sorted_dp[valid[-n:]], sorted_idx[valid[-n:]], gd.x[sorted_idx[valid[-n:]]], gd.y[sorted_idx[valid[-n:]]]])
    with open('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/large_elev.txt', 'w') as file:  
        file.write('x y z id\n')
        np.savetxt(file, np.c_[gd.x[sorted_idx[valid[-n:]]], gd.y[sorted_idx[valid[-n:]]], sorted_dp[valid[-n:]], sorted_idx[valid[-n:]]])
    with open('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/small_elev.txt', 'w') as file:  
        file.write('x y z id\n')
        np.savetxt(file, np.c_[gd.x[sorted_idx[valid[:n]]], gd.y[sorted_idx[valid[:n]]], sorted_dp[valid[:n]], sorted_idx[valid[:n]]])


def quality_check_hgrid(gd, epsg=4326, outdir='./'):
    '''
    Check several types of grid issues that may crash the model:
    '''

    if epsg != 26918:
        gd.proj(prj0=f'epsg:{epsg}', prj1='epsg:26918')

    # bad quads
    bp_name = f'{outdir}/bad_quad.bp'
    gd.check_quads(angle_min=60,angle_max=120,fname=bp_name)
    bad_quad_bp = Bpfile(filename=bp_name)
    if bad_quad_bp.n_nodes > 0:
        new_gr3_name = f'hgrid_split_quads_{epsg}.gr3'
        gd.split_quads(angle_min=60,angle_max=120,fname=f'{outdir}/{new_gr3_name}')
        gd = schism_grid(f'{outdir}/{new_gr3_name}')
        print(f'{bad_quad_bp.n_nodes} bad quads splitted and the updated hgrid is saved as {new_gr3_name}')

    # small and skew elements
    gd.compute_all()
    sorted_area = np.sort(gd.area)
    sorted_idx = np.argsort(gd.area)

    small_ele = sorted_area < 5.0
    print(f'\n{sum(small_ele)} small (< 5.0 m2) elements:')
    print(np.c_[sorted_area[small_ele], sorted_idx[small_ele]+1, gd.xctr[sorted_idx[small_ele]], gd.yctr[sorted_idx[small_ele]]])
    if sum(small_ele) > 0:
        Bpfile(xyz_array=np.c_[gd.xctr[sorted_idx[small_ele]], gd.yctr[sorted_idx[small_ele]], sorted_idx[small_ele]+1]).writer(f'{outdir}/small_ele.bp')

    # skew elements
    bp_name = f'{outdir}/skew_ele.bp'
    skew_ele = gd.check_skew_elems(angle_min=0.5, fmt=1, fname=bp_name) 
    SMS_MAP(detached_nodes=np.c_[gd.xctr[skew_ele], gd.yctr[skew_ele], gd.yctr[skew_ele]*0]).writer(f'{outdir}/skew_ele.map')

    small_ele = np.argwhere(gd.area < 5.0).reshape(-1, )
    SMS_MAP(detached_nodes=np.c_[gd.xctr[small_ele], gd.yctr[small_ele], gd.yctr[small_ele]*0]).writer(f'{outdir}/small_ele.map')

    invalid = np.unique(np.array([*skew_ele, *small_ele]))
    if len(invalid) > 0:
        invalid_elnode = gd.elnode[invalid].reshape(-1, 1)
        invalid_elnode = np.unique(invalid_elnode)
        invalid_elnode = invalid_elnode[invalid_elnode>=0]
        invalid_neighbors = np.unique(gd.ine[invalid_elnode].reshape(-1, 1))
        invalid_neighbors = invalid_neighbors[invalid_neighbors>=0]
        
        SMS_MAP(detached_nodes=np.c_[gd.xctr[invalid_neighbors], gd.yctr[invalid_neighbors], gd.yctr[invalid_neighbors]*0]).writer(f'{outdir}/invalid_element_relax.map')
        return gd

def pre_proc_hgrid(hgrid_name=''):

    dirname = os.path.dirname(hgrid_name)
    file_basename = os.path.basename(hgrid_name)
    file_extension = pathlib.Path(hgrid_name).suffix

    if file_extension == '.2dm':
        gd = sms2grd(hgrid_name)
        gd.write_hgrid(f'{dirname}/hgrid.utm.26918')
    elif file_extension == '.gr3' or file_extension == '.ll':
        gd = schism_grid(hgrid_name)
    else:
        raise Exception('Extension unknown')

    quality_check_hgrid(gd, epsg=26918, outdir=dirname)

    gd.proj(prj0='epsg:26918', prj1='epsg:4326')
    gd.save(f'{dirname}/hgrid.ll')

    # load bathymetry
    # print('submit bathymetry loading script to queue')
    # os.chdir(dirname)
    # os.system('./pload_depth.py')
    # print('done loading bathymetry')
    # time.sleep(10000000)
    pass

def tweak_depths(hgrid_name=''):

    dirname = os.path.dirname(hgrid_name)
    file_basename = os.path.basename(hgrid_name)
    file_extension = pathlib.Path(hgrid_name).suffix

    os.system(f'mv {hgrid_name} {dirname}/hgrid.ll')
    gd = schism_grid(f'{dirname}/hgrid.ll')
    gd_DEM_loaded = copy.deepcopy(gd)

    print('loading levee heights from National Levee Database')
    gd = set_levee_profile(hgrid=gd, wdir=dirname, levee_info_dir=f'{dirname}/Levee_info/')
    print('loading additional tweaks on levee heights')
    gd = set_additional_dp_v11_91(gd_ll=gd, gd_dem=gd_DEM_loaded, wdir=dirname, levee_info_dir=f'{dirname}/Levee_info/Additional_Polygons/')

    print('outputing hgrid.ll')
    gd.save(f'{dirname}/hgrid.ll')

def gen_hgrid_formats(hgrid_name='', gd:schism_grid=None):
    if gd is None:
        gd = read_schism_hgrid_cached(hgrid_name, overwrite_cache=True)
    else:
        hgrid_name = gd.source_file

    gd.lon, gd.lat = gd.x, gd.y

    dirname = os.path.dirname(hgrid_name)
    file_basename = os.path.basename(hgrid_name)
    file_extension = pathlib.Path(hgrid_name).suffix

    if os.popen(f'grep "open" {dirname}/hgrid.ll').read() == '':
        os.system(f'cat {dirname}/bnd >> {dirname}/hgrid.ll')

    print('outputing hgrid.cpp')
    gd_cpp = copy.deepcopy(gd)
    gd_cpp.x, gd_cpp.y = lonlat2cpp(lon=gd.x, lat=gd.y, lon0=-77.07, lat0=24.0)
    gd_cpp.save(f'{dirname}/hgrid.cpp.gr3')
    os.system(f'mv {dirname}/hgrid.cpp.gr3 {dirname}/hgrid.cpp')
    if os.popen(f'grep "open" {dirname}/hgrid.cpp').read() == '':
        os.system(f'cat {dirname}/bnd >> {dirname}/hgrid.cpp')
    gd.cpp_x, gd.cpp_y = gd_cpp.x, gd_cpp.y

    print('outputing hgrid in UTM')
    gd.proj(prj0='epsg:4326', prj1='epsg:26918')
    gd.write_hgrid(f'{dirname}/hgrid.utm.26918.gr3')
    gd.utm_x, gd.utm_y = gd.x, gd.y
    if os.popen(f'grep "open" {dirname}/hgrid.utm.26918.gr3').read() == '':
        os.system(f'cat {dirname}/bnd >> {dirname}/hgrid.utm.26918.gr3')
    print('outputing *.2dm')
    grd2sms(gd, f'{dirname}/hgrid.utm.2dm')

    print('outputing hgrid.102008.gr3')
    gd.x, gd.y = gd.lon, gd.lat
    gd.proj(prj0='epsg:4326', prj1='esri:102008')
    gd.write_hgrid(f'{dirname}/hgrid.102008.gr3')
    gd.x_102008, gd.y_102008 = gd.x, gd.y
    if os.popen(f'grep "open" {dirname}/hgrid.102008.gr3').read() == '':
        os.system(f'cat {dirname}/bnd >> {dirname}/hgrid.102008.gr3')

    print('saving *.pkl, which has x, y of all needed projections')
    with open(f'{dirname}/hgrids.pkl', 'wb') as file:
        pickle.dump(gd, file)

    with open(f'{dirname}/hgrids.pkl', 'rb') as file:
       gd_test = pickle.load(file)
    print('finish generating hgrids in different projections/formats')
    pass

if __name__ == "__main__":
    # pre_proc_hgrid('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/v14.35_relaxed2.2dm')
    # tweak_depths('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/hgrid.ll.new')
    gen_hgrid_formats('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23l/Hgrid/hgrid.ll')
    
    '''
    gd_fname = '/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23p11/hgrid.ll'
    gd_cache_fname = os.path.splitext(gd_fname)[0] + '.pkl'
    if os.path.exists(gd_cache_fname):
        gd = schism_grid(gd_cache_fname)
    else:
        gd = schism_grid(gd_fname)
        gd.save(gd_cache_fname)
    
    quality_check_hgrid(gd)
    '''
    pass
