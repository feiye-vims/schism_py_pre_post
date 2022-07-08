# %%
import enum
from matplotlib.image import AxesImage
from matplotlib.pyplot import ylabel
from pylib import loadz, proj_pts
from schism_py_pre_post.Grid.SMS import get_all_points_from_shp, SMS_ARC, SMS_MAP, curvature
import numpy as np
import copy
import os
import math
from osgeo import gdal
from dataclasses import dataclass
import pathlib
import pickle


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

def improve_thalwegs(S, x, y, xt_left, yt_left, xt_right, yt_right, search_steps=200):
    xts = np.linspace(xt_left, xt_right, search_steps, axis=1)
    yts = np.linspace(yt_left, yt_right, search_steps, axis=1)

    jj, ii = Sidx(S, xts[:], yts[:])
    elevs = S.elev[ii, jj]
    real_thalweg_idx = np.argmin(elevs, axis=1)

    x_real = xts[np.array(range(len(xts))), real_thalweg_idx]
    y_real = yts[np.array(range(len(yts))), real_thalweg_idx]
    z_real = elevs[np.array(range(len(yts))), real_thalweg_idx]
    
    return [x_real, y_real]

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
    line_cplx = np.squeeze(line.view(np.complex128))
    dist = np.absolute(line_cplx[1:] - line_cplx[:-1])

    # return np.r_[0.0, np.cumsum(dist)]
    return np.r_[dist[0], dist]

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

def smooth_bank(line, xs, ys, xs_other_side, ys_other_side, ang_diff_shres=np.pi/2.4, nmax=100):
    n = 0
    while n < nmax:
        angle_diffs = np.r_[0.0, get_angle_diffs(xs, ys), 0.0]
        sharp_turns = np.argwhere(abs(angle_diffs) > ang_diff_shres)[:, 0]
        if not sharp_turns.size:
            break
        else:
            # Step 1: plan to move the sharp turn point to the average coordinates of the two adjacent points
            xs_moved = np.array((xs[sharp_turns-1] + xs[sharp_turns+1]) / 2)
            ys_moved = np.array((ys[sharp_turns-1] + ys[sharp_turns+1]) / 2)

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

    SMS_MAP(arcs=range_arcs).writer(filename=f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_test_range.map')

    return x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, bank2bank_width

def redistribute_arc(line, channel_width):
    # along-river distance, for redistribution
    dist_along_thalweg = get_dist_increment(line)

    curv = curvature(line)
    for i in [2, 3]:
        for j in range(i):
            curv[j::i] = np.maximum(curv[j::i], curvature(line[j::i]))

    R = 1.0/(curv+1e-10)
    river_resolution = np.minimum(0.3 * R, 5 * channel_width)
    river_resolution = np.minimum(np.maximum(10, river_resolution), 300)

    retained_points = np.ones((dist_along_thalweg.shape), dtype=bool)
    idx = 1
    original_seg_length = dist_along_thalweg[1]
    while idx < len(dist_along_thalweg)-1:
        ideal_seg_length = river_resolution[idx]  # move ideal arc forward
        if original_seg_length < ideal_seg_length:
            retained_points[idx] = False  # remove points
            original_seg_length += dist_along_thalweg[idx]
        else:
            original_seg_length = dist_along_thalweg[idx]
        idx += 1  # move original arc forward
    # last point should be retained
    retained_points[-1] = True
            
    return line[retained_points, :]

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
        
def assemble_arcs(arcs, line, blast_radius, thalweg_id, i_check_valid=False):
    valid_idx = np.ones(len(line[:, 0]), dtype=bool)
    if i_check_valid:
        for k in [0, -1]:
            valid_idx *= (line[:, 0] - line[k, 0])**2 + (line[:, 1] - line[k, 1])**2 > blast_radius[k]**2
        if sum(valid_idx) < 1:
            print(f'warning: thalweg {thalweg_id+1} has less than 1 points after bombing, neglecting ...')
            return arcs  # no changes

    arcs.append(SMS_ARC(points=np.c_[line[valid_idx], line[valid_idx]]))
    return arcs

if __name__ == "__main__":
    MapUnit2METER = 1 #00000

    # tif_fname = r'/sciclone/data10/wangzg/DEM/npz/sc_ll_7.npz'  # directory of DEM data

    # tif_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/CUDEM_TX_merged_LA_utm15.tif'
    # thalweg_shp_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_1_0_fix_region1.shp'

    tif_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_1.0_Fix_region2_utm15N.tif'
    thalweg_shp_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_1_0_fix_region2_test1.shp'

    river_threshold = np.array([15, 400]) / MapUnit2METER
    nudge_ratio = np.array((0.3, 2.0))  # ratio between nudging distance to mean half-channel-width
    
    # %%
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

    # read points from thalweg shapefile
    xyz, l2g, curv = get_all_points_from_shp(thalweg_shp_fname)
    thalwegs = []
    thalwegs_curv = []
    thalweg_endpoints = np.empty((0, 2), dtype=float)
    for i, idx in enumerate(l2g):
        print(f'Arc {i+1} of {len(l2g)}')
        thalwegs.append(xyz[idx, :])
        thalwegs_curv.append(curv[idx])
        thalweg_endpoints = np.r_[thalweg_endpoints, np.reshape(xyz[idx][0, :], (1, 2))]
        thalweg_endpoints = np.r_[thalweg_endpoints, np.reshape(xyz[idx][-1, :], (1, 2))]

    # Dry run
    thalweg_endpoints_width = np.empty((0), dtype=float)
    thalweg_widths = []
    valid_thalwegs = []
    for i, [line, curv] in enumerate(zip(thalwegs, thalwegs_curv)):
        print(f'Dry run: Arc {i+1} of {len(l2g)}')
        x_banks_left, _, x_banks_right, _, _, width = get_two_banks(line)
        thalweg_widths.append(width)
        if width is None:
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, 0.0]
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, 0.0]
        else:
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, width[0]]
            thalweg_endpoints_width = np.r_[thalweg_endpoints_width, width[1]]
        if len(line[:, 0]) < 2:
            print(f"warning: thalweg {i+1} only has one point, neglecting ...")
            valid_thalwegs.append(False)
            continue

        if x_banks_left is None or x_banks_right is None:
            print(f"warning: thalweg {i+1} out of DEM coverage, neglecting ...")
            valid_thalwegs.append(False)
            continue

        valid_thalwegs.append(True)
    # End Dry run
    
    bank_arcs = []
    bank_arcs_all = []
    redistributed_arcs = []
    # improved_thalweg_arcs = []
    for i, [line, curv, width, valid_thalweg] in enumerate(zip(thalwegs, thalwegs_curv, thalweg_widths, valid_thalwegs)):
        print(f'Arc {i+1} of {len(l2g)}')

        if not valid_thalweg:
            print(f"marked as invalid in dry run, skipping ...")
            continue

        # redistribute vertices
        line = redistribute_arc(line, width)
        redistributed_arcs.append(SMS_ARC(points=np.c_[line[:, 0], line[:, 1]]))

        if len(line[:, 0]) < 2:
            print(f"warning: thalweg {i+1} only has one point after redistribution, neglecting ...")
            continue

        # re-make banks based on redistributed thalweg
        x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, width = get_two_banks(line)

        # touch-ups on the two banks
        if x_banks_left is None or x_banks_right is None:
            print(f'warning: cannot find banks for thalweg {i+1} after redistribution, neglecting ... ')
            continue

        # nudge banks
        x_banks_left, y_banks_left = nudge_bank(line, perp+np.pi, x_banks_left, y_banks_left, dist=nudge_ratio*0.5*np.mean(width))
        x_banks_right, y_banks_right = nudge_bank(line, perp, x_banks_right, y_banks_right, dist=nudge_ratio*0.5*np.mean(width))

        # smooth banks
        line, x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp = smooth_bank(line, x_banks_left, y_banks_left, x_banks_right, y_banks_right)
        line, x_banks_right, y_banks_right, x_banks_left, y_banks_left, perp = smooth_bank(line, x_banks_right, y_banks_right, x_banks_left, y_banks_left)

        # make inner arcs between two banks
        x_inner_arcs = np.linspace(x_banks_left, x_banks_right, 4)
        y_inner_arcs = np.linspace(y_banks_left, y_banks_right, 4)

        # remove points at thalweg intersectios with a radius of channel width
        blast_radius = np.array([0.0, 0.0])
        for k in [0, -1]:
            dist = thalweg_endpoints[:, :] - line[k, :]
            neighbor_thalwegs_endpoints = np.argwhere(dist[:, 0]**2 + dist[:, 1]**2 < 1e-6**2)
            if len(neighbor_thalwegs_endpoints) > 1:
                blast_radius[k] = 0.4 * np.mean(thalweg_endpoints_width[neighbor_thalwegs_endpoints])
            else:  # no intersections
                blast_radius[k] = 0.0

        out_arcs = [
            np.c_[x_banks_left, y_banks_left],
            np.c_[x_banks_right, y_banks_right],
        ]
        for x_inner_arc, y_inner_arc in zip (x_inner_arcs, y_inner_arcs):
            out_arcs.append(np.c_[x_inner_arc, y_inner_arc])

        for out_arc in out_arcs:
            bank_arcs = assemble_arcs(bank_arcs, out_arc, blast_radius, i, i_check_valid=True)
            bank_arcs_all = assemble_arcs(bank_arcs_all, out_arc, blast_radius, i, i_check_valid=False)

    # write SMS maps
    SMS_MAP(arcs=bank_arcs).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/bank.map')
    SMS_MAP(arcs=bank_arcs_all).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/bank_all.map')
    SMS_MAP(arcs=redistributed_arcs).writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/redist_thalweg.map')

    # %%
    pass
