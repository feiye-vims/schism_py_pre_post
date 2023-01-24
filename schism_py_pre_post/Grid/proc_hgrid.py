# %%
from codecs import ignore_errors
from logging import raiseExceptions
from fasteners import InterProcessLock
from pylib import schism_grid, sms2grd, grd2sms
import shutil
import os
import numpy as np
from pyrsistent import get_in
from schism_py_pre_post.Grid.Bpfile import Bpfile
from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached, get_inp, propogate_nd, get_bnd_nd_cached
from schism_py_pre_post.Grid.SMS import SMS_MAP, lonlat2cpp, cpp2lonlat
from schism_py_pre_post.Shared_modules.set_levee_height import set_constant_levee_height
from schism_py_pre_post.Shared_modules.set_levee_profile import set_levee_profile
from schism_py_pre_post.Shared_modules.set_additional_dp import set_additional_dp_v11_91
from schism_py_pre_post.Shared_modules.set_feeder_dp import set_feeder_dp
import pathlib
import copy
import pickle
import lloyd
from scipy import spatial
import subprocess



def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[1]

def compute_area(gd, ie):
    fp=gd.elnode[ie,-1]<0;
    x1=gd.x[gd.elnode[ie,0]]; y1=gd.y[gd.elnode[ie,0]];
    x2=gd.x[gd.elnode[ie,1]]; y2=gd.y[gd.elnode[ie,1]];
    x3=gd.x[gd.elnode[ie,2]]; y3=gd.y[gd.elnode[ie,2]];
    x4=gd.x[gd.elnode[ie,3]]; y4=gd.y[gd.elnode[ie,3]]; x4[fp]=x1[fp]; y4[fp]=y1[fp]
    area=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)+(x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2
    return area

def spring(gd):
    return gd

def reduce_bad_elements(
    gd=None, fixed_eles_id=np.array([], dtype=int), fixed_points_id=np.array([], dtype=int),
    area_threshold=150, skewness_threshold=10, nmax=8, output_fname=None
):

    n = 0
    while True:
        n += 1
        # make a copy of the grid in case the procedure below fails
        gd0 = copy.deepcopy(gd)

        gd.compute_bnd()
        gd.compute_ic3()

        if n==1:
            fixed_points_id = np.sort(np.unique(np.r_[fixed_points_id, gd.bndinfo.ip]))
        else:
            fixed_points_id = gd.bndinfo.ip

        gd.compute_nne()
        fixed_eles_id = np.unique(gd.ine[fixed_points_id].flatten())
        fixed_eles_id = fixed_eles_id[fixed_eles_id>=0]
        i_fixed_eles = np.zeros((gd.ne, ), dtype=bool)
        i_fixed_eles[fixed_eles_id] = True

        reduce_type = 2 + n % 2  #alternating two methods

        gd.compute_area()
        gd.compute_nne()

        is_small_triangle = (gd.area < area_threshold) * (gd.i34 == 3)
        small_triangle_id = np.argwhere(is_small_triangle)

        small_elnode_x = gd.x[gd.elnode[is_small_triangle, :3]]
        small_elnode_y = gd.y[gd.elnode[is_small_triangle, :3]]
        small_elnode_xy = small_elnode_x[:, :3] + 1j*small_elnode_y[:, :3]
        
        small_elnode_flatten = gd.elnode[is_small_triangle, :3].flatten()
        _, idx, counts = np.unique(small_elnode_xy.flatten(), return_counts=True, return_index=True)
        degenerate_count = np.zeros((gd.np, ), dtype=int)
        degenerate_count[small_elnode_flatten[idx]] = counts

        small_triangle_dl = np.zeros((sum(is_small_triangle), 3), dtype=float)
        for i in range(3):
            j = (i + 1) % 3
            small_triangle_dl[:, i] = abs(small_elnode_xy[:, i] - small_elnode_xy[:, j])

        small_triangle_dl_sorted = np.sort(small_triangle_dl, axis=1)
        aspect_ratios =  -np.log(abs(1.0 - small_triangle_dl_sorted[:, 2] / np.sum(small_triangle_dl_sorted[:, :2], axis=1)))

        gd_xy = gd.x + 1j*gd.y

        i_ele_shorted = np.zeros((gd.ne, ), dtype=bool)

        if reduce_type == 3:
            # ------------3 point reduce --------------------------------
            small_triangle_reduce = aspect_ratios < skewness_threshold
        else:
            small_triangle_reduce = aspect_ratios >= skewness_threshold
            gd.compute_side(fmt=2)
        small_triangle_reduce_id = small_triangle_id[small_triangle_reduce]

        if not small_triangle_reduce.any():
            print(f'Done fixing type ({reduce_type}) bad elements in {n-1} iterations')
            break
        elif i_fixed_eles[small_triangle_reduce_id].all():
            print(f'Done fixing type ({reduce_type}) bad elements in {n-1} iterations, except for fixed eles: {small_triangle_reduce_id}')
            break
        elif n > 1 and n <= nmax:
            leftover = np.zeros((gd.ne,), dtype=bool)
            leftover[small_triangle_id] = True
            leftover[fixed_eles_id] = False
            print(f'Failed to fix {sum(leftover)} eles: \n{np.argwhere(leftover).flatten()[:50]+1}\n after {n-1} iterations')
            if n == nmax:
                break

        print(f'\n>>>>Element reduction Iteration {n}: trying to fix {len(small_triangle_id)-sum(i_fixed_eles[small_triangle_id])} eles; reduce_type: {reduce_type}\n')
        print(f'min element area = {np.sort(gd.area)[:20]}')

        i_ele_reduce = np.zeros((gd.ne, ), dtype=bool)
        i_ele_reduce[small_triangle_reduce_id] = True

        idx = np.argsort(gd.area[small_triangle_id].flatten())
        small_triangle_id_sorted = small_triangle_id[idx].flatten()
        aspect_ratio_sorted = aspect_ratios[idx].flatten()

        for i, [ie, aspect_ratio] in enumerate(zip(small_triangle_id_sorted, aspect_ratio_sorted)):
            if i_ele_shorted[ie] or i_fixed_eles[ie]:
                continue
            if reduce_type == 3 and aspect_ratio >= skewness_threshold:
                continue
            elif reduce_type == 2 and aspect_ratio < skewness_threshold:
                continue

            iee = ie
            for j in range(2):
                iee = np.unique(gd.ine[gd.elnode[iee]].flatten())
                iee = iee[iee>=0]
            i_ele_shorted[iee] = True

            if reduce_type == 3:
                this_sorted_nodes = np.sort(gd.elnode[ie, :3].flatten())
                this_degenerate_count = degenerate_count[this_sorted_nodes]
                i_degenerate_node = np.zeros((3, ), dtype=bool)
                iee = [x for x in iee if not (x in gd.ic3[ie] or x==ie)]
            elif reduce_type == 2:  # find the shortest side as the degenerate side
                elsides = gd.elside[ie]
                dl = gd.distj[elsides[elsides>=0]]
                degenerate_side = np.argmin(dl)
                this_sorted_nodes = np.sort(gd.isidenode[elsides[degenerate_side], :])
                this_degenerate_count = degenerate_count[this_sorted_nodes]
                i_degenerate_node = np.zeros((2, ), dtype=bool)
                iee = [x for x in iee if not (x == gd.ic3[ie, degenerate_side] or x==ie)]

            i_degenerate_node[np.argmax(this_degenerate_count)] = True

            gd_xy0 = copy.deepcopy(gd_xy)
            gd_xy[this_sorted_nodes[~i_degenerate_node]] = \
                gd.x[this_sorted_nodes[i_degenerate_node]] + 1j*gd.y[this_sorted_nodes[i_degenerate_node]]  # reduce the last two points to the first point
            gd.x = np.real(gd_xy)
            gd.y = np.imag(gd_xy)
            area = compute_area(gd, iee)
            if min(area) < 1e-12:
                i_ele_shorted[iee] = False
                i_ele_shorted[ie] = True
                gd_xy = copy.deepcopy(gd_xy0)


        # ----- update basic grid info for output ------
        # update new node sequence
        gd_xy_unique, inv = np.unique(gd_xy, return_inverse=True)
        # update nodes
        gd.x = np.real(gd_xy_unique)
        gd.y = np.imag(gd_xy_unique)
        gd.dp = np.zeros((len(gd.x),), dtype=float)

        # update elements
        inv = np.r_[inv, -2, -2]  # pad two elements because elnode uses -2 as null
        # map node ids in ele map to new node sequence
        gd.elnode = inv[gd.elnode]

        # identify degenerate trianlges with duplicated node ids (new sequence)
        duplicate_nd_eles = np.zeros((gd.elnode.shape[0], ), dtype=bool)
        duplicate_nd_eles += (gd.elnode[:, 0] == gd.elnode[:, 1])
        duplicate_nd_eles += (gd.elnode[:, 0] == gd.elnode[:, 2])
        duplicate_nd_eles += (gd.elnode[:, 0] == gd.elnode[:, 3])
        duplicate_nd_eles += (gd.elnode[:, 1] == gd.elnode[:, 2])
        duplicate_nd_eles += (gd.elnode[:, 1] == gd.elnode[:, 3])
        duplicate_nd_eles += (gd.elnode[:, 2] == gd.elnode[:, 3])

        # remove degenerate ele
        gd.elnode = gd.elnode[~duplicate_nd_eles]
        gd.i34 = gd.i34[~duplicate_nd_eles]

        i_fixed_eles = i_fixed_eles[~duplicate_nd_eles]

        gd.ne = gd.elnode.shape[0]
        gd.np = gd.x.shape[0]

        gd.compute_area()
        if min(gd.area) < 0.0:
            print(f'found negative elements: {np.argwhere(gd.area<0.0)+1}')
            grd2sms(gd, 'failed.2dm')
            print(f'final element areas: {np.sort(gd.area)}')
            print(f'reverting to previous grid ...')
            gd = copy.deepcopy(gd0)
        else:
            tmp_filename = f'{os.path.dirname(output_fname)}/{pathlib.Path(output_fname).stem}.fix{n}.2dm'
            grd2sms(gd, tmp_filename)
            gd = sms2grd(tmp_filename)
            pass

    # end of while loop

    if output_fname is None:
        gd.save('tmp.gr3')
        gd = schism_grid('tmp.gr3')
        os.remove('tmp.gr3')
    else:
        if pathlib.Path(output_fname).suffix == '.2dm':
            grd2sms(gd, output_fname)
            gd = sms2grd(output_fname)
        elif pathlib.Path(output_fname).suffix in ['.gr3', 'll']:
            gd.save(output_fname)
            gd = schism_grid(output_fname)
        else:
            raiseExceptions('suffix of output grid file name not supported.')

    return gd

def lloyd_relax(gd, target_points):
    # import matplotlib.pyplot as plt
    # plt.triplot(new_xy[:,0], new_xy[:,1], tri.simplices)
    # plt.plot(new_xy[:,0], new_xy[:,1], 'o')
    # plt.show()

    return gd
# %%  
def grid_element_relax(gd, target_points=None, niter=3, ntier=0, max_dist=50, min_area_allowed=1e-3, wdir=None, output_fname=None):
    if not os.path.exists(wdir):
        raise Exception(f'wdir: {wdir} not found')

    # set fixed points (which won't be moved)
    if target_points is not None:
        if target_points.dtype == bool:
            pass
        elif target_points.dtype == int:
            tmp = np.zeros((gd.np, ), dtype=bool)
            tmp[target_points] = True
            target_points = tmp
        else:
            raise Exception('unsupported fixe_points type')
    else:  # all points can be moved if target_points are not specified
        target_points = np.ones((gd.np, ), dtype=bool)

    # fix boundary points too
    interior_points = np.ones((gd.np, ), dtype=bool)
    gd.compute_bnd()
    interior_points[gd.bndinfo.ip] = False

    gd.write_hgrid(f'{wdir}/hgrid_spring_input.gr3', fmt=1)

    # find inp
    if not hasattr(gd, 'ine'):
        gd.compute_nne()

    i_relax = target_points * interior_points
    relax_points = np.argwhere(i_relax).flatten()

    # # Lloyd
    #     field = lloyd.Field(np.c_[gd.x[nd_nei], gd.y[nd_nei]])
    #     field.relax()
    #     new_xy = field.get_points()

    # Optional: multiple tiers (limited by max_dist)
    inp2 = get_inp(gd, ntiers=ntier).reshape(gd.np, -1)
    gd_xy = np.r_[gd.x + 1j*gd.y, 1e10+1j*1e10]
    nd_ids = np.array(range(gd.np))
    inp_distance = abs(gd_xy[inp2] - np.repeat(gd_xy[nd_ids].reshape(-1, 1), inp2.shape[1], axis=1))
    inp_shortrange = copy.deepcopy(inp2)
    inp_shortrange[inp_distance>max_dist] = -1

    relax_points = np.unique(inp_shortrange[target_points].flatten())
    relax_points = relax_points[relax_points>=0]
    i_relax = np.zeros((gd.np, ), dtype=bool)
    i_relax[relax_points] = True
    i_relax = i_relax * interior_points
    relax_points = np.argwhere(i_relax)
    
    gd_fixed = copy.deepcopy(gd)
    gd_fixed.dp[:] = 0
    gd_fixed.dp[~i_relax] = 1
    gd_fixed.save(f'{wdir}/fixed.gr3')
    
    # springing
    p = subprocess.Popen(f'{wdir}/grid_spring', cwd=wdir, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    p.stdin.write(f'{niter}\n{min_area_allowed}\n1e6\n1\n'.encode()) #expects a bytes type object
    p.communicate()[0]
    p.stdin.close()

    # check negative element and revert points to original locations
    # n = 0
    # while True:
    #     n += 1
    #     for i in range(niter):
    #         print (f'iteration {i+1} of {niter}')
    #         gd_x_padded = copy.deepcopy(np.r_[gd.x, np.nan])
    #         gd_y_padded = copy.deepcopy(np.r_[gd.y, np.nan])

    #         gd.x = np.nanmean(gd.x[inp[i_relax]], axis=1)
    #         gd.y = np.nanmean(gd.y[inp[i_relax]], axis=1)

    #         if n == 10:
    #             grd2sms(gd, f'{os.path.dirname(gd.source_file)}/hgrid_failed.2dm')
    #             raise Exception(f'Negative area cannot be fixed during relaxing: {np.argwhere(negative_ele)+1}')

    #         gd.compute_area()
    #         negative_ele = gd.area < 1e-10
    #         if any(negative_ele):
    #             if n == 1:  # revert all nodes of the negative elements
    #                 revert_nodes = np.unique(gd.elnode[negative_ele])
    #                 revert_nodes = revert_nodes[revert_nodes>=0]
    #             else:  # revert all nodes of the negative elements, as well as their neighboring nodes, tiers expanding with each try
    #                 for i in range(n):
    #                     revert_nodes = np.unique(inp[revert_nodes])
    #                     revert_nodes = revert_nodes[revert_nodes>=0]
    #             gd.x[revert_nodes] = gd_x0[revert_nodes]
    #             gd.y[revert_nodes] = gd_y0[revert_nodes]
    #         else:
    #             break

    gd = schism_grid(f'{wdir}/hgrid.spring.gr3')
    grd2sms(gd, f'{wdir}/hgrid.spring.2dm')

    if output_fname is not None:
        if pathlib.Path(output_fname).suffix == '.2dm':
            shutil.move(f'{wdir}/hgrid.spring.2dm', output_fname)
        elif pathlib.Path(output_fname).suffix in ['.gr3', 'll']:
            shutil.move(f'{wdir}/hgrid.spring.gr3', output_fname)
        else:
            raiseExceptions('suffix of output grid file name not supported.')

    gd.compute_area()
    print(f'Sorted area after springing: {np.sort(gd.area)[:50]}')

    return gd
    
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


def quality_check_hgrid(gd, outdir='./', small_ele_limit=5.0, skew_ele_minangle=0.5):
    '''
    Check several types of grid issues that may crash the model:
    If input grid is lon/lat, reproject it so that the projected unit is meter
    '''

    print('\n-----------------quality check>')
    # small and skew elements
    gd.compute_all()
    sorted_area = np.sort(gd.area)
    sorted_idx = np.argsort(gd.area)

    # small elements
    i_small_ele = sorted_area < small_ele_limit
    n_small_ele = sum(i_small_ele)
    print(f'\n{n_small_ele} small (< {small_ele_limit} m2) elements:')
    small_ele = np.argwhere(i_small_ele).flatten()
    if n_small_ele > 0:
        print(np.c_[sorted_area[small_ele], sorted_idx[small_ele]+1, gd.xctr[sorted_idx[small_ele]], gd.yctr[sorted_idx[small_ele]]][:min(10, n_small_ele)])
        # Bpfile(xyz_array=np.c_[gd.xctr[sorted_idx[small_ele]], gd.yctr[sorted_idx[small_ele]], sorted_idx[small_ele]+1]).writer(f'{outdir}/small_ele.bp')
        # small_ele = np.argwhere(gd.area < small_ele_limit).reshape(-1, )
        # SMS_MAP(detached_nodes=np.c_[gd.xctr[small_ele], gd.yctr[small_ele], gd.yctr[small_ele]*0]).writer(f'{outdir}/small_ele.map')

    # skew elements
    bp_name = f'{outdir}/skew_ele.bp'
    skew_ele = gd.check_skew_elems(angle_min=skew_ele_minangle, fmt=1, fname=bp_name) 
    print(f'\n{len(skew_ele)} skew (min angle < {skew_ele_minangle})')
    if len(skew_ele) > 0:
        print(skew_ele[:min(10, len(skew_ele))])
        # SMS_MAP(detached_nodes=np.c_[gd.xctr[skew_ele], gd.yctr[skew_ele], gd.yctr[skew_ele]*0]).writer(f'{outdir}/skew_ele.map')

    invalid = np.unique(np.array([*skew_ele, *small_ele]))
    if len(invalid) > 0:
        invalid_elnode = gd.elnode[invalid].reshape(-1, 1)
        invalid_elnode = np.unique(invalid_elnode)
        invalid_elnode = invalid_elnode[invalid_elnode>=0]
        invalid_neighbors = np.unique(gd.ine[invalid_elnode].reshape(-1, 1))
        invalid_neighbors = invalid_neighbors[invalid_neighbors>=0]
        i_invalid_nodes = np.zeros((gd.np, ), dtype=bool)
        i_invalid_nodes[invalid_elnode] = True  # set invalid nodes to be "not fixed", i.e., can be tweaked.
        
        # SMS_MAP(detached_nodes=np.c_[gd.xctr[invalid_neighbors], gd.yctr[invalid_neighbors], gd.yctr[invalid_neighbors]*0]).writer(f'{outdir}/invalid_element_relax.map')
    else:
        invalid_elnode = None
        invalid_neighbors = None
        i_invalid_nodes = None


    return {'hgrid': gd, 'invalid_nodes': invalid_elnode, 'invalid_elements': invalid_neighbors, 'i_invalid_nodes': i_invalid_nodes}

def pre_proc_hgrid(hgrid_name='', prj='esri:102008', load_bathy=False, nmax=4):
    '''
    Fix small and skew elements and bad quads
    prj: needs to specify hgrid's projection (the unit must be in meters)
    nmax: maximum number of rounds of fixing, most fixable elements can be fixed with nmax=4,
          nmax>4 ususally doesn't make additional improvements
    '''
    prj_name = prj.replace(':', '_')

    gd, dir_info = read_schism_hgrid_cached(hgrid_name, return_source_dir=True)
    dirname = dir_info['dir']
    file_basename = dir_info['basename']
    file_extension = dir_info['extension']

    # Fix invalid elements
    n_fix = 0
    while True:
        n_fix += 1

        grid_quality = quality_check_hgrid(gd, outdir=dirname)
        i_target_nodes = grid_quality['i_invalid_nodes']

        if i_target_nodes is None:  # all targets fixed
            print('\n -------------------------Done fixing invalid elements --------------------------------------------')
            break
        elif n_fix > nmax:  # maximum iteration reached, exit with leftovers
            print(' --------------------------------------------Done fixing invalid elements,')
            print(f"but failed at the following elements: {grid_quality['invalid_elements']}")
            print(f"and nodes: {grid_quality['invalid_nodes']}")
            break
        else:  # fix targets

            print(f'\n----------------Fixing invalid elements, Round {n_fix}--------------------')

            # split bad quads
            print('\n ------------------- Splitting bad quads >')
            bp_name = f'{dirname}/bad_quad.bp'
            gd.check_quads(angle_min=60,angle_max=120,fname=bp_name)
            bad_quad_bp = Bpfile(filename=bp_name)
            if bad_quad_bp.n_nodes > 0:
                new_gr3_name = f"hgrid_split_quads.gr3"
                gd.split_quads(angle_min=60,angle_max=120,fname=f'{dirname}/{new_gr3_name}')
                gd = schism_grid(f'{dirname}/{new_gr3_name}')
                print(f'{bad_quad_bp.n_nodes} bad quads split and the updated hgrid is saved as {new_gr3_name}')

                # quality check again since gd is updated
                i_target_nodes = quality_check_hgrid(gd, outdir=dirname)['i_invalid_nodes']
                if i_target_nodes is None: continue

            if n_fix == 1:
                # include intersection relax points only at Round 1
                inter_relax_pts = SMS_MAP(filename=f'{dirname}/intersection_relax.map').detached_nodes
                gd.xy = gd.x + 1j* gd.y
                gd.proj(prj0=prj, prj1='epsg:4326')
                _, inter_relax_nd = spatial.cKDTree(np.c_[gd.x, gd.y]).query(inter_relax_pts[:, :2])
                gd.x, gd.y = np.real(gd.xy), np.imag(gd.xy)

                i_target_nodes[inter_relax_nd] = True

            print('\n ------------------- Reducing small/skew elements>')
            print(f'Number of target nodes: {sum(i_target_nodes)}')
            target_nodes_expand, i_target_nodes_expand = propogate_nd(gd, i_target_nodes, ntiers=3)
            gd = reduce_bad_elements(
                gd=gd, fixed_points_id=np.argwhere(~i_target_nodes_expand).flatten(),
                area_threshold=80,
                output_fname=f'{dirname}/{file_basename}_fix_bad_eles_round_{n_fix}.2dm'
            )

            # quality check again since gd is updated 
            i_target_nodes = quality_check_hgrid(gd, outdir=dirname)['i_invalid_nodes']
            if i_target_nodes is None: continue

            if n_fix == 1: 
                # re-find intersection nodes since gd is updated
                gd.xy = gd.x + 1j* gd.y
                gd.proj(prj0=prj, prj1='epsg:4326')
                _, inter_relax_nd = spatial.cKDTree(np.c_[gd.x, gd.y]).query(inter_relax_pts[:, :2])
                gd.x, gd.y = np.real(gd.xy), np.imag(gd.xy)

                i_target_nodes[inter_relax_nd] = True


            print('\n ------------------- Relaxing remaining small/skew elements>')
            gd = grid_element_relax(
                gd=gd, target_points=i_target_nodes, niter=min(10, n_fix), ntier=2, max_dist=50, wdir=dirname,
                output_fname=f'{dirname}/{file_basename}_relax_round_{n_fix}.2dm'
            )

            print(f'\n****************Done fixing invalid elements, Round {n_fix}*********************\n')

    grd2sms(gd, f'{dirname}/hgrid.2dm')
    gd.proj(prj0=prj, prj1='epsg:4326')
    gd.save(f'{dirname}/hgrid.ll')

    # load bathymetry
    # if load_bathy:
    #     os.chdir(dirname)
    #     os.system('python ./pload_depth.py')
    #     if not os.exist(f'{dirname}/hgrid.ll.new'):
    #         raise Exception('failed to load bathymetry')
    #     else:
    #         print('done loading bathymetry')
            
    pass

def tweak_depths(hgrid_name=''):

    dirname = os.path.dirname(hgrid_name)
    file_basename = os.path.basename(hgrid_name)
    file_extension = pathlib.Path(hgrid_name).suffix

    gd_ll_original = schism_grid(f'{dirname}/hgrid.ll')
    

    os.system(f'mv {hgrid_name} {dirname}/hgrid.ll')
    gd = schism_grid(f'{dirname}/hgrid.ll')
    gd.x, gd.y = gd_ll_original.x, gd_ll_original.y
    gd_DEM_loaded = copy.deepcopy(gd)

    print('set default levee heights (-9 m)')
    gd = set_constant_levee_height(gd=gd, wdir=dirname)

    print('loading levee heights from National Levee Database')
    gd = set_levee_profile(gd=gd, wdir=dirname)

    print('loading additional tweaks on levee heights')
    gd = set_additional_dp_v11_91(gd_ll=gd, gd_dem=gd_DEM_loaded, wdir=dirname)

    print('outputing hgrid.ll')
    gd.save(f'{dirname}/hgrid.ll')

    # set feeder dp
    gd = set_feeder_dp(
        feeder_info_dir='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/',
        new_grid_dir=dirname
    )

    os.system(f'mv {dirname}/hgrid.ll {dirname}/hgrid.ll_before_feeder_dp')
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
    # Sample usage

    # wdir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/v14.42_post_proc2/'
    wdir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/Relax_test4/'

    # Step 1: check grid quality.
    # pre_proc_hgrid(f'{wdir}/v14.42.gr3')
    pre_proc_hgrid(f'{wdir}/Relax_test4.2dm')

    # Step 1.5
    # load DEM using pload

    # Step 2
    tweak_depths(f'{wdir}/hgrid.ll.new')  # renamed to hgrid.ll after this step

    # Step 3
    gen_hgrid_formats(f'{wdir}/hgrid.ll')  # put bnd in the wdir before this step
    pass