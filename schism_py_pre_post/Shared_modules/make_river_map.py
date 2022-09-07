# %%
import matplotlib.pyplot as plt
from pylib import loadz, proj_pts
from schism_py_pre_post.Grid.SMS import get_all_points_from_shp, SMS_ARC, SMS_MAP, curvature
import numpy as np
import copy
from copy import deepcopy
import os
import math
from osgeo import gdal
from dataclasses import dataclass
import pathlib
import pickle
from scipy import interpolate


@dataclass
class dem_data():
    x: np.ndarray
    y: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    elev: np.ndarray
    dx: float
    dy: float

@dataclass
class Bomb():
    x: np.ndarray
    y: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    parent_thalwegs : list
    blast_radius: np.ndarray


# %%
def getAngle(a, b, c):
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return ang + 360 if ang < 0 else ang

# %%
def Sidx(S, x, y):
    ''' return nearest index (i, j) in DEM mesh for point (x, y) '''
    dSx = S.lon[1] - S.lon[0]
    dSy = S.lat[1] - S.lat[0]
    i = (np.round((x - S.lon[0]) / dSx)).astype(int)
    j = (np.round((y - S.lat[0]) / dSy)).astype(int)
    return [i, j]

def improve_thalwegs(S, line, search_length, perp):
    x = line[:, 0]
    y = line[:, 1]

    xt_right = line[:, 0] + search_length * np.cos(perp)
    yt_right = line[:, 1] + search_length * np.sin(perp)
    xt_left = line[:, 0] + search_length * np.cos(perp + np.pi)
    yt_left = line[:, 1] + search_length * np.sin(perp + np.pi)

    __search_steps = int(np.max(search_length/ds))

    j, i = Sidx(S, x, y)
    I_lr = np.c_[i, i] * 0
    J_lr = copy.deepcopy(I_lr)
    thalweg_idx = np.ones(xt_right.shape) * 9999
    for k, [xt, yt] in enumerate([[xt_left, yt_left], [xt_right, yt_right]]):
        xts = np.linspace(x, xt, __search_steps, axis=1)
        yts = np.linspace(y, yt, __search_steps, axis=1)

        jj, ii = Sidx(S, xts[:], yts[:])
        try:
            elevs = S.elev[ii, jj]
            low = np.argpartition(elevs, min(10, elevs.shape[1]-1), axis=1)
            thalweg_idx = np.median(low[:, :10], axis=1).astype(int)
            # thalweg_idx = np.argmin(elevs, axis=1)

            I_lr[:, k] = ii[range(0, len(x)), thalweg_idx]
            J_lr[:, k] = jj[range(0, len(x)), thalweg_idx]
        except IndexError:
            return None, None

    lr = S.elev[I_lr[:,0], J_lr[:,0]] > S.elev[I_lr[:,1], J_lr[:,1]]
    j = J_lr[range(0, len(x)), lr.astype(int)]
    i = I_lr[range(0, len(x)), lr.astype(int)]
    
    x_real = S.lon[j]
    y_real = S.lat[i]

    return np.c_[x_real, y_real]

# %%
def get_bank(S, x, y, eta, xt, yt, search_steps=100, search_tolerance=5):
    '''Get a bank on one side of the thalweg (x, y)'''
    # search_steps_tile = np.repeat(np.arange(search_steps).reshape(1, -1), len(x), axis=0)  # expanded to the search area

    # form a search area between thalweg and search limit
    xts = np.linspace(x, xt, search_steps, axis=1)
    yts = np.linspace(y, yt, search_steps, axis=1)

    j, i = Sidx(S, x, y)
    elevs_thalweg = np.ones(S.elev[i, j].shape) * eta  # elev on thalweg
    elevs_stream = np.tile(elevs_thalweg.reshape(-1, 1), (1, search_steps))  # expanded to the search area

    jj, ii = Sidx(S, xts[:], yts[:])
    try:
        elevs = S.elev[ii, jj]
    except IndexError:
        return None, None

    R = (elevs - elevs_stream)  # closeness to target depth: 0-1
    bank_idx = np.argmax(R>0, axis=1)
    
    invalid = bank_idx == 0
    bank_idx[invalid] = np.argmin(abs(R[invalid, :]), axis=1)

    # R_sort_idx = np.argsort(R)
    # bank_idx = np.min(R_sort_idx[:, :min(search_steps, search_tolerance)], axis=1)

    x_banks = S.lon[jj[range(0, len(x)), bank_idx]]
    y_banks = S.lat[ii[range(0, len(x)), bank_idx]]

    return x_banks, y_banks

def get_dist_increment(line):
    line_copy = copy.deepcopy(line)
    line_cplx = np.squeeze(line_copy.view(np.complex128))
    dist = np.absolute(line_cplx[1:] - line_cplx[:-1])

    # return np.r_[0.0, np.cumsum(dist)]
    return dist

def get_angle_diffs(xs, ys):
    line = np.c_[xs, ys]
    line_cplx = np.squeeze(line.view(np.complex128))
    angles = np.angle(np.diff(line_cplx))
    angle_diff0 = np.diff(angles)
    angle_diff = np.diff(angles)
    angle_diff[angle_diff0 > np.pi] -= 2 * np.pi
    angle_diff[angle_diff0 < -np.pi] += 2 * np.pi

    return angle_diff

def ccw(A,B,C):
    # A is a point with the coordinate x=A[0], y=A[1]
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def smooth_thalweg(line, ang_diff_shres=np.pi/2.4, nmax=100, smooth_coef=0.2):
    xs = line[:, 0]
    ys = line[:, 1]

    n = 0
    while n < nmax:
        angle_diffs = np.r_[0.0, get_angle_diffs(xs, ys), 0.0]
        sharp_turns = np.argwhere(abs(angle_diffs) > ang_diff_shres)[:, 0]
        if not sharp_turns.size:
            break
        else:
            # Step 1: plan to move the sharp turn point to a new point, based on
            # the average coordinates of the two adjacent points and the original location
            line[sharp_turns, 0] = np.array((xs[sharp_turns-1] + xs[sharp_turns+1]) / 2) * smooth_coef + xs[sharp_turns] * (1-smooth_coef)
            line[sharp_turns, 1] = np.array((ys[sharp_turns-1] + ys[sharp_turns+1]) / 2) * smooth_coef + ys[sharp_turns] * (1-smooth_coef)

            n += 1

    if n == nmax:
        print(f'warning: smooth_thalweg did not converge in {n} steps\n')
    else:
        print(f'smooth_thalweg converged in {n} steps\n')

    perp = get_perpendicular_angle(line)
        
    return line, perp

def river_quality(xs, ys, idx, ang_diff_shres=np.pi/2.4):
    # identify channels that are too ambiguous and hopeless
    if sum(idx) < 2:
        return False
        
    for [x, y] in zip(xs, ys):
        angle_diffs = abs(np.r_[0.0, get_angle_diffs(x[idx], y[idx]), 0.0])
        if np.sum(angle_diffs)/len(x[idx]) > 0.75:
            print(f'discarding arc based on the number of sharp turns\n')
            return False
    
    return True

def smooth_bank(line, xs, ys, xs_other_side, ys_other_side, ang_diff_shres=np.pi/2.4, nmax=100, smooth_coef=0.2):
    n = 0
    while n < nmax:
        angle_diffs = np.r_[0.0, get_angle_diffs(xs, ys), 0.0]
        sharp_turns = np.argwhere(abs(angle_diffs) > ang_diff_shres)[:, 0]

        if len(sharp_turns)==0:
            break
        else:
            # Step 1: plan to move the sharp turn point to a new point, based on
            # the average coordinates of the two adjacent points and the original location
            xs_moved = np.array((xs[sharp_turns-1] + xs[sharp_turns+1]) / 2) * smooth_coef + xs[sharp_turns] * (1-smooth_coef)
            ys_moved = np.array((ys[sharp_turns-1] + ys[sharp_turns+1]) / 2) * smooth_coef + ys[sharp_turns] * (1-smooth_coef)

            # Step 2: decide if the planned move is too far (i.e., across the thalweg)
            insert_turns_idx = []
            for i, [x_moved, y_moved, sharp_turn] in enumerate(zip(xs_moved, ys_moved, sharp_turns)):
                invalid_move = intersect([xs[sharp_turn], ys[sharp_turn]], [x_moved, y_moved], line[sharp_turn, :], line[sharp_turn-1, :]) + \
                               intersect([xs[sharp_turn], ys[sharp_turn]], [x_moved, y_moved], line[sharp_turn, :], line[sharp_turn+1, :])
                if invalid_move:
                    # prepare to insert 2 more points around the sharp turn
                    insert_turns_idx.append(i) 
                else:
                    xs[sharp_turns] = xs_moved
                    ys[sharp_turns] = ys_moved
            
            if len(insert_turns_idx) > 0:
                # replace the point at sharp bend with two more points adjacent to it
                idx = sharp_turns[insert_turns_idx]

                tmp = np.c_[line, xs, ys, xs_other_side, ys_other_side]
                tmp_original = copy.deepcopy(tmp)
                tmp[idx, :] = (tmp_original[idx, :] + tmp_original[idx+1, :])/2
                tmp = np.insert(tmp, idx, (tmp_original[idx, :] + tmp_original[idx-1, :])/2, axis=0)

                line = tmp[:, :2]
                xs = tmp[:, 2]
                ys = tmp[:, 3]
                xs_other_side = tmp[:, 4]
                ys_other_side = tmp[:, 5]

            n += 1

    if n == nmax:
        print(f'warning: smooth_bank did not converge in {n} steps\n')
    else:
        print(f'smooth_bank converged in {n} steps\n')

    perp = get_perpendicular_angle(line)
        
    return line, xs, ys, xs_other_side, ys_other_side, perp

def nudge_bank(line, perp, xs, ys, dist=np.array([35, 500])):
    ds = ((line[:, 0] - xs)**2 + (line[:, 1] - ys)**2)**0.5

    idx = ds < dist[0]
    xs[idx] = line[idx, 0] + dist[0] * np.cos(perp[idx])
    ys[idx] = line[idx, 1] + dist[0] * np.sin(perp[idx])

    idx = ds > dist[1]
    xs[idx] = line[idx, 0] + dist[1] * np.cos(perp[idx])
    ys[idx] = line[idx, 1] + dist[1] * np.sin(perp[idx])

    return xs, ys

def Tif2XYZ(tif_fname=None):
    cache_name = tif_fname + '.pkl'
    if os.path.exists(cache_name):
        with open(cache_name, 'rb') as f:
            S = pickle.load(f)
    else:
        ds = gdal.Open(tif_fname, gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)

        width = ds.RasterXSize
        height = ds.RasterYSize

        gt = ds.GetGeoTransform()
        TL_x, TL_y = gt[0], gt[3]

        #showing a 2D image of the topo
        # plt.imshow(elevation, cmap='gist_earth',extent=[minX, maxX, minY, maxY])
        # plt.show()

        z = band.ReadAsArray()

        dx = gt[1]
        dy = gt[5]
        if gt[2] != 0 or gt[4] != 0:
            raise Exception()

        x_idx = np.array(range(width))
        y_idx = np.array(range(height))
        xp = dx * x_idx + TL_x + dx/2
        yp = dy * y_idx + TL_y + dy/2

        S = dem_data(xp, yp, xp, yp, z, dx, dy)
        with open(cache_name, 'wb') as f:
            pickle.dump(S, f, protocol=pickle.HIGHEST_PROTOCOL)

    return S

def get_perpendicular_angle(line):
    line_cplx = np.squeeze(line.copy().view(np.complex128))
    angles = np.angle(np.diff(line_cplx))
    angle_diff0 = np.diff(angles)
    angle_diff = np.diff(angles)
    angle_diff[angle_diff0 > np.pi] -= 2 * np.pi
    angle_diff[angle_diff0 < -np.pi] += 2 * np.pi
    perp = angles[:-1] + angle_diff / 2 - np.pi / 2
    perp = np.r_[angles[0] - np.pi / 2, perp, angles[-1] - np.pi / 2]

    return perp

def set_eta(x, y):
    # thalweg_eta = np.maximum(0.0, (y - 3313760.0))/(3367300.0 - 3313760.0) * 1.2
    # thalweg_eta = np.ones(y.shape) * 0.5

    y0 = [0, 23388517, 3461404, 9e9]
    eta0 = [0, 0, 0, 0]

    eta = np.interp(y, y0, eta0)

    return eta

class river_arc():
    def __init__(self, line=None):
        self.thalweg = line
        self.perp = get_perpendicular_angle(line)
        self.endpoints = np.array([line[0, :], line[-1, :]])

class river_seg():
    def __init__(self, thalweg: river_arc):
        self.thalweg = thalweg
    
    def make_inner_arcs(self, narcs=2):
        inner_arcs = np.linspace(self.left_bank, self.right_bank, narcs)

def get_two_banks(line):
    range_arcs = []
    # find perpendicular direction along thalweg at each point
    perp = get_perpendicular_angle(line)

    # find search area for a thalweg, consisting of two lines on each side
    xt_right = line[:, 0] + search_length * np.cos(perp)
    yt_right = line[:, 1] + search_length * np.sin(perp)
    xt_left = line[:, 0] + search_length * np.cos(perp + np.pi)
    yt_left = line[:, 1] + search_length * np.sin(perp + np.pi)

    # Diagnostic: save search area as SMS arcs
    range_arcs += [SMS_ARC(points=np.c_[xt_left, yt_left]), SMS_ARC(points=np.c_[xt_right, yt_right])]

    # set water level at each point along the thalweg, based on observation, simulation, estimation, etc.
    thalweg_eta = set_eta(line[:, 0], line[:, 1])

    # find two banks
    x_banks_left, y_banks_left = get_bank(S, line[:, 0], line[:, 1], thalweg_eta, xt_left, yt_left, search_steps=search_steps)
    x_banks_right, y_banks_right = get_bank(S, line[:, 0], line[:, 1], thalweg_eta, xt_right, yt_right, search_steps=search_steps)

    # get attributes of the initial banks
    # average width, for deciding nudging distance
    if x_banks_left is None or x_banks_right is None:
        print('warning: failed to find banks ... ')
        return None, None, None, None, None, None

    bank2bank_width = ( (x_banks_left - x_banks_right)**2 + (y_banks_left - y_banks_right)**2 ) **0.5

    SMS_MAP(arcs=range_arcs).writer(filename=f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/bank_range.map')

    return x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, bank2bank_width

def moving_average(a, n=10, self_weights=0):
    if len(a) <= n:
        ret2 = a * 0.0 + np.mean(a)
        return ret2
    else:
        ret = np.cumsum(a, axis=0, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ret[n-1:] = ret[n-1:] / n

        # re-align time series
        ret1 = ret * 0.0
        m = int(np.floor(n/2))
        ret1[m:-m] = ret[2*m:]

        # fill the first and last few records
        ret1[:m] = ret1[m]
        ret1[-m:] = ret1[-m-1]

        # put more weights on self
        ret2 = (ret1 + self_weights * a) / (1 + self_weights)

        return ret2

def redistribute_arc(line, line_smooth, channel_width, smooth_option=1, R_coef=0.4, width_coef=4.0, reso_thres=[3, 300]):
    # along-river distance, for redistribution
    dist_along_thalweg = get_dist_increment(line)

    # # smooth curvature to get rid of small-scale zig-zags
    if smooth_option == -1:  # existing line_smooth
        curv = curvature(line_smooth)
    elif smooth_option == 0:  # no smoothing
        line_smooth = line
        curv = curvature(line)
    elif smooth_option == 1:  # Option 1: moving average
        line_smooth = moving_average(line, n=30, self_weights=2)
        curv = curvature(line_smooth)
    elif smooth_option == 2:  # Option 2: spline (slow and doesn't really work because original points are preserved)
        smooth_factor = 4
        #create spline function
        f, u = interpolate.splprep([line[:, 0], line[:, 1]], s=10, per=0)
        #create interpolated lists of points
        uint = np.interp(np.arange(0, len(u)-1+1/smooth_factor, 1/smooth_factor), np.arange(0, len(u)), u)
        xint, yint = interpolate.splev(uint, f)
        line_smooth = np.c_[xint, yint]
        curv_sp = curvature(line_smooth)
        curv = curv_sp[::smooth_factor]
    '''
    plt.scatter(line_smooth[:, 0], line_smooth[:, 1])
    plt.scatter(line[:, 0], line[:, 1], s=0.3)
    plt.show()
    '''

    # use different interval along the line to calculate curvature
    if smooth_option != 2:
        for i in [1, 2]:  # more combinations -> more conservative, i.e., larger curvature
            for j in range(i):
                curv[j::i] = np.maximum(curv[j::i], curvature(line_smooth[j::i]))

    R = 1.0/(curv+1e-10)
    river_resolution = np.minimum(R_coef * R, width_coef * channel_width)
    river_resolution = np.minimum(np.maximum(reso_thres[0], river_resolution), reso_thres[1])
    river_resolution_seg = (river_resolution[:-1]+river_resolution[1:])/2  # resolution between two points

    retained_points = np.ones((line.shape[0]), dtype=bool)
    idx = 0
    this_seg_length = dist_along_thalweg[0]  # dist between pt0 and pt1
    while idx < len(dist_along_thalweg)-1:
        if this_seg_length < river_resolution_seg[idx]:  # resolution of the seg between pt0 and pt1
            retained_points[idx+1] = False  # remove point
            this_seg_length += dist_along_thalweg[idx+1]
        else:
            this_seg_length = dist_along_thalweg[idx+1]
        idx += 1  # move original arc forward
    # last point should be retained
    retained_points[-1] = True
            
    return line[retained_points, :], line_smooth, river_resolution, retained_points

def snap_vertices(line, thalweg_resolution):
    dist_along_thalweg = get_dist_increment(line)

    idx = 0
    original_seg_length = dist_along_thalweg[0]
    while idx < len(dist_along_thalweg)-1:
        if original_seg_length < thalweg_resolution[idx]:
            line[idx+1, :] = line[idx, :]  # snap current point to the previous one
            original_seg_length += dist_along_thalweg[idx+1]
        else:
            original_seg_length = dist_along_thalweg[idx+1]
        idx += 1  # move original arc forward

    return line

def nearest_neighbour(points_a, points_b):
    from scipy import spatial
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[0], tree.query(points_a)[1]

def get_thalweg_neighbors(thalwegs, thalweg_endpoints):
    # rounding_digits = 12
    # thalweg_endpoints_unique = np.unique(np.array(thalweg_endpoints).round(decimals=rounding_digits), axis=0)
    # thalweg_bombs = np.zeros((len(thalwegs), 2), dtype=bool)

    thalweg_neighbors = [None] * len(thalwegs) * 2
    for i, thalweg in enumerate(thalwegs):
        dist = thalweg_endpoints[:, :] - thalweg[0, :]
        same_points = dist[:, 0]**2 + dist[:, 1]**2 < 1e-10**2
        if sum(same_points) > 1:  # intersection point
            thalweg_neighbors[2*i] = np.argwhere(same_points)[0]
    
    return thalweg_neighbors
        
def bomb_line(line, blast_radius, thalweg_id, i_check_valid=False):
    valid_idx = np.ones(len(line[:, 0]), dtype=bool)
    if i_check_valid:
        for k in [0, -1]:
            valid_idx *= (line[:, 0] - line[k, 0])**2 + (line[:, 1] - line[k, 1])**2 > blast_radius[k]**2
        if sum(valid_idx) < 1:
            print(f'warning: thalweg {thalweg_id+1} has less than 1 points after bombing, neglecting ...')
            return valid_idx  # no changes

    return valid_idx

if __name__ == "__main__":

    # ------------------------- basic inputs --------------------------- 
    MapUnit2METER = 1 #00000

    # tif_fname = r'/sciclone/data10/wangzg/DEM/npz/sc_ll_7.npz'  # directory of DEM data

    # tif_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/CUDEM_TX_merged_LA_utm15.tif'
    # thalweg_shp_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_1_0_fix_region1.shp'

    # tif_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_1.0_Fix_region2_utm15N.tif'
    # thalweg_shp_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_1_0_fix_region2_test1.shp'

    # tif_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_dem_merged_utm17N.tif'
    # thalweg_shp_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_riverstreams_cleaned_utm17N.shp'
    # thalweg_smooth_shp_fname = None  # '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_riverstreams_cleaned_corrected_utm17N.shp'

    tif_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_dem_merged_utm17N.tif'
    thalweg_shp_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/v4.shp'
    thalweg_smooth_shp_fname = None  # '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_riverstreams_cleaned_corrected_utm17N.shp'

    river_threshold = np.array([15, 400]) / MapUnit2METER
    nudge_ratio = np.array((0.3, 2.0))  # ratio between nudging distance to mean half-channel-width

    blast_radius_scale = 0.4  # coef controlling the blast radius at intersections
    intersect_res_scale  = 0.4  # coef controlling the resolution of the paved mesh at intersections
    # ------------------------- end basic inputs --------------------------- 
    
    # %%
    # ------------------------- read DEM --------------------------- 
    if pathlib.Path(tif_fname).suffix == ".npz":
        S = loadz(tif_fname)
    elif pathlib.Path(tif_fname).suffix == ".tif" : 
        S = Tif2XYZ(tif_fname=tif_fname)
    else:
        raise Exception("Unknown DEM format.")
    print(f'DEM box: {min(S.lon)}, {min(S.lat)}, {max(S.lon)}, {max(S.lat)}')

    dx = S.lon[1] - S.lon[0]
    dy = S.lat[1] - S.lat[0]
    ds = (abs(dx) + abs(dy)) / 2

    search_length = river_threshold[1] * 1.1
    search_steps = int(river_threshold[1] / ds)

    starndard_watershed_resolution = 400.0  # meters
    nrow_arcs = 4  # the channel is resolved by "nrow_arcs" rows of elements
    # ------------------------- read thalweg --------------------------- 
    xyz, l2g, curv = get_all_points_from_shp(thalweg_shp_fname)
    if thalweg_smooth_shp_fname is not None:
        xyz_s, l2g_s, curv_s = get_all_points_from_shp(thalweg_smooth_shp_fname)

    thalwegs = []
    thalwegs_smooth = []
    thalwegs_curv = []
    thalweg_endpoints = np.empty((0, 2), dtype=float)
    for i, idx in enumerate(l2g):
        # print(f'Arc {i+1} of {len(l2g)}')
        thalwegs.append(xyz[idx, :])
        thalwegs_curv.append(curv[idx])
        thalweg_endpoints = np.r_[thalweg_endpoints, np.reshape(xyz[idx][0, :], (1, 2))]
        thalweg_endpoints = np.r_[thalweg_endpoints, np.reshape(xyz[idx][-1, :], (1, 2))]
        if thalweg_smooth_shp_fname is not None:
            thalwegs_smooth.append(xyz_s[idx, :])
        else:
            thalwegs_smooth.append(None)

    # ------------------------- Dry run ---------------------------
    print('Dry run')
    thalweg_endpoints_width = np.empty((0), dtype=float)
    thalweg_widths = []
    valid_thalwegs = []
    original_banks = []
    for i, [line, curv] in enumerate(zip(thalwegs, thalwegs_curv)):
        # print(f'Dry run: Arc {i+1} of {len(l2g)}')
        x_banks_left, y_banks_left, x_banks_right, y_banks_right, _, width = get_two_banks(line)
        thalweg_widths.append(width)
        if width is None:
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, 0.0]
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, 0.0]
        else:
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, width[0]]
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, width[-1]]

        if len(line[:, 0]) < 2:
            print(f"warning: thalweg {i+1} only has one point, neglecting ...")
            valid_thalwegs.append(False)
            continue

        if x_banks_left is None or x_banks_right is None:
            print(f"warning: thalweg {i+1} out of DEM coverage, neglecting ...")
            valid_thalwegs.append(False)
            continue

        original_banks.append(SMS_ARC(points=np.c_[x_banks_left, y_banks_left]))
        original_banks.append(SMS_ARC(points=np.c_[x_banks_right, y_banks_right]))
        valid_thalwegs.append(True)

    # End Dry run: found valid river segments; record channel width
    
    # ------------------------- Wet run --------------------------- 
    # initialize some lists and array to hold the arc information
    bank_arcs = np.empty((len(thalwegs), 2), dtype=object) # left bank and right bank for each thalweg
    bank_arcs_raw = deepcopy(bank_arcs)
    inner_arcs = np.empty((len(thalwegs), nrow_arcs), dtype=object)
    cc_arcs = deepcopy(bank_arcs)  # [, 0] is head, [, 1] is tail
    smoothed_thalwegs = [None] * len(thalwegs)
    redistributed_thalwegs = [None] * len(thalwegs)
    corrected_thalwegs = [None] * len(thalwegs)
    final_thalwegs = [None] * len(thalwegs)
    intersection_res_scatters = []
    thalwegs_neighbors = deepcopy(bank_arcs)  # [, 0] is head, [, 1] is tail
    real_bank_width = np.zeros((len(thalwegs), 2), dtype=float)  # [, 0] is head, [, 1] is tail

    # enumerate each thalweg
    for i, [thalweg, curv, width, valid_thalweg, thalweg_smooth] in enumerate(zip(thalwegs, thalwegs_curv, thalweg_widths, valid_thalwegs, thalwegs_smooth)):
        print(f'Wet run: Arc {i+1} of {len(l2g)}')

        if not valid_thalweg:
            print(f"marked as invalid in dry run, skipping ...")
            continue

        # Redistribute thalwegs vertices
        thalweg, thalweg_smooth, reso, retained_idx = redistribute_arc(thalweg, thalweg_smooth, width, smooth_option=1)
        smoothed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg_smooth[:, 0], thalweg_smooth[:, 1]])
        redistributed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]])

        if len(thalweg[:, 0]) < 2:
            print(f"warning: thalweg {i+1} only has one point after redistribution, neglecting ...")
            continue

        # re-make banks based on redistributed thalweg
        x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, width = get_two_banks(thalweg)

        # correct thalwegs
        width_moving_avg = moving_average(width, n=10)
        thalweg = improve_thalwegs(S, thalweg, width_moving_avg*0.5, perp)
        corrected_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]])

        # Redistribute thalwegs vertices
        thalweg, thalweg_smooth, reso, retained_idx = redistribute_arc(thalweg, thalweg_smooth[retained_idx], width, smooth_option=1)
        smoothed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg_smooth[:, 0], thalweg_smooth[:, 1]])
        redistributed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]])

        # Smooth thalweg
        thalweg, perp = smooth_thalweg(thalweg, ang_diff_shres=np.pi/2.4)
        final_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]])

        # re-make banks based on corrected thalweg
        x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, width = get_two_banks(thalweg)

        # touch-ups on the two banks
        if x_banks_left is None or x_banks_right is None:
            print(f'warning: cannot find banks for thalweg {i+1} after redistribution, neglecting ... ')
            continue

        # nudge banks
        x_banks_left, y_banks_left = nudge_bank(thalweg, perp+np.pi, x_banks_left, y_banks_left, dist=nudge_ratio*0.5*np.mean(width))
        x_banks_right, y_banks_right = nudge_bank(thalweg, perp, x_banks_right, y_banks_right, dist=nudge_ratio*0.5*np.mean(width))

        # smooth banks
        thalweg, x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp = smooth_bank(thalweg, x_banks_left, y_banks_left, x_banks_right, y_banks_right)
        if thalweg is None:
            continue
        thalweg, x_banks_right, y_banks_right, x_banks_left, y_banks_left, perp = smooth_bank(thalweg, x_banks_right, y_banks_right, x_banks_left, y_banks_left)
        if thalweg is None:
            continue

        # update width
        width = ((x_banks_left-x_banks_right)**2 + (y_banks_left-y_banks_right)**2)**0.5
        
        # get actual resolution along redistributed/smoothed thalweg
        thalweg_resolution = get_dist_increment(thalweg)

        # make inner arcs between two banks
        x_inner_arcs = np.linspace(x_banks_left, x_banks_right, nrow_arcs)
        y_inner_arcs = np.linspace(y_banks_left, y_banks_right, nrow_arcs)

        # determine blast radius based on mean channel width at an intersection
        blast_radius = np.array([0.0, 0.0])
        for k in [0, -1]:  # head and tail
            dist = thalweg_endpoints[:, :] - thalweg[k, :]
            neighbor_thalwegs_endpoints = np.argwhere(dist[:, 0]**2 + dist[:, 1]**2 < 200**2)
            thalwegs_neighbors[i, k] = neighbor_thalwegs_endpoints
            if len(neighbor_thalwegs_endpoints) > 1:
                blast_radius[k] = blast_radius_scale * np.mean(thalweg_endpoints_width[neighbor_thalwegs_endpoints])
            else:  # no intersections
                blast_radius[k] = 0.0

        # assemble banks
        valid_points = bomb_line(np.c_[x_banks_left, y_banks_left], blast_radius, i, i_check_valid=True) * \
                       bomb_line(np.c_[x_banks_right, y_banks_right], blast_radius, i, i_check_valid=True)
        for k, line in enumerate([np.c_[x_banks_left, y_banks_left], np.c_[x_banks_right, y_banks_right]]):
            bank_arcs_raw[i, k] = SMS_ARC(points=np.c_[line[:, 0], line[:, 1]])
            if sum(valid_points) > 0:
                bank_arcs[i, k] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1]])

        if sum(valid_points) > 0:
            for j in [0, -1]:
                real_bank_width[i, j] = ((x_banks_left[valid_points][j]-x_banks_right[valid_points][j])**2 + (y_banks_left[valid_points][j]-y_banks_right[valid_points][j])**2)**0.5

        # quality check river arcs
        if river_quality(x_inner_arcs, y_inner_arcs, valid_points):
            # assemble inner arcs
            for k, [x_inner_arc, y_inner_arc] in enumerate(zip(x_inner_arcs, y_inner_arcs)):
                line = np.c_[x_inner_arc, y_inner_arc]
                if sum(valid_points) > 0:
                    line = snap_vertices(line, width * 0.3)  # optional: thalweg_resolution*0.75
                    inner_arcs[i, k] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1]])

            # assemble cross-channel arcs
            if sum(valid_points) > 0:
                for j in [0, -1]:
                    cc_arcs[i, j] = SMS_ARC(points=np.c_[x_inner_arcs[:, valid_points][:, j], y_inner_arcs[:, valid_points][:, j]])

    # assemble intersectional resolution scatters
    for i, thalweg_neighbors in enumerate(thalwegs_neighbors):
        for j, neibs in enumerate(thalweg_neighbors):  # head and tail
            if neibs is not None and len(neibs) > 1:
                # intersect_ring is the effective refining region, consisting of bank endpoints at the intersection of multiple thalwegs
                intersect_ring = np.zeros((0, 3), dtype=float)
                for j, nei in enumerate(neibs):
                    id = int(nei/2)
                    i_tail_head = int(nei%2)
                    if bank_arcs[id, 0] is not None:
                        intersect_ring = np.r_[intersect_ring, np.c_[bank_arcs[id, 0].nodes[i_tail_head, :2].reshape(1, -1), intersect_res_scale * real_bank_width[id, i_tail_head]]]
                    if bank_arcs[id, 1] is not None:
                        intersect_ring = np.r_[intersect_ring, np.c_[bank_arcs[id, 1].nodes[i_tail_head, :2].reshape(1, -1), intersect_res_scale * real_bank_width[id, i_tail_head]]]

                # make a scatter set consisting of a center point (where further refining is optional) and a outer ring (where resolution reverts to normal value)
                if len(intersect_ring) > 0:
                    center_point = np.mean(intersect_ring, axis=0).reshape(1, -1)

                    diff = intersect_ring - center_point
                    mean_radius = np.mean((diff[:, 0]**2 + diff[:, 1]**2)**0.5)

                    # outter ring is nudged to 120% based on the intersection ring
                    outter_ring = (intersect_ring - center_point) * 1.2 + center_point
                    # nudge the outter ring points to make the outter ring more or less a circle, avoiding outter_ring points inside intersection ring
                    diff = outter_ring - center_point
                    radius = (diff[:, 0]**2 + diff[:, 1]**2)**0.5
                    stretch = (np.maximum(1.2, radius/mean_radius)/(radius/mean_radius)).reshape(-1, 1)
                    outter_ring = np.tile(stretch, (1,3)) * (outter_ring - center_point) + center_point

                    outter_ring[:, -1] = starndard_watershed_resolution

                    # assemble all rings
                    all = np.r_[intersect_ring, center_point, outter_ring]
                    intersection_res_scatters.append(all)

    intersect_res = None
    if intersection_res_scatters:  # len > 0
        intersect_res = np.concatenate(intersection_res_scatters, axis=0)
        valid = intersect_res[:, -1] > 0
        intersect_res = intersect_res[valid, :]
        
    # assemble river arcs
    # river_arcs = np.r_[bank_arcs.reshape((-1, 1)), inner_arcs.reshape((-1, 1))]
    # End wet run

    # ------------------------- write SMS maps --------------------------- 
    SMS_MAP(arcs=bank_arcs.reshape((-1, 1))).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/bank.map')
    SMS_MAP(arcs=cc_arcs.reshape((-1, 1))).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/cc_arcs.map')
    SMS_MAP(arcs=bank_arcs_raw.reshape((-1, 1))).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/bank_raw.map')
    SMS_MAP(arcs=inner_arcs.reshape((-1, 1))).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/river.map')
    SMS_MAP(arcs=smoothed_thalwegs).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/smoothed_thalweg.map')
    SMS_MAP(arcs=redistributed_thalwegs).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/redist_thalweg.map')
    SMS_MAP(arcs=corrected_thalwegs).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/corrected_thalweg.map')
    SMS_MAP(arcs=final_thalwegs).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/final_thalweg.map')
    if intersect_res is not None:
        np.savetxt('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/intersection_res.xyz', intersect_res)

    # %%
    pass
