# %%
from matplotlib.image import AxesImage
from pylib import loadz
from schism_py_pre_post.Grid.SMS import get_all_points_from_shp, SMS_ARC, SMS_MAP
import numpy as np
import copy
import os
import math
from osgeo import gdal
from dataclasses import dataclass
import pathlib


@dataclass
class dem_data():
     x: np.ndarray
     y: np.ndarray
     lon: np.ndarray
     lat: np.ndarray
     elev: np.ndarray
     dx: float
     dy: float


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
    elevs = S.elev[ii, jj]

    R = abs(elevs - elevs_stream)  # closeness to target depth: 0-1
    R_sort_idx = np.argsort(R)
    bank_idx = np.min(R_sort_idx[:, :min(search_steps, search_tolerance)], axis=1)

    x_banks = S.lon[jj[range(0, len(x)), bank_idx]]
    y_banks = S.lat[ii[range(0, len(x)), bank_idx]]

    return x_banks, y_banks

def get_angle_diffs(xs, ys):
    line = np.c_[xs, ys]
    line_cplx = np.squeeze(line.view(np.complex128))
    angles = np.angle(np.diff(line_cplx))
    angle_diff0 = np.diff(angles)
    angle_diff = np.diff(angles)
    angle_diff[angle_diff0 > np.pi] -= 2 * np.pi
    angle_diff[angle_diff0 < -np.pi] += 2 * np.pi

    return angle_diff

def smooth_bank(xs, ys, ang_diff_shres=np.pi/2.5, nmax=100):
    n = 0
    while n < nmax:
        angle_diffs = np.r_[0.0, get_angle_diffs(xs, ys), 0.0]
        sharp_turns = np.squeeze(np.argwhere(abs(angle_diffs) > ang_diff_shres))
        if not sharp_turns.size:
            break
        else:
            # xs[sharp_turns] = xs[sharp_turns-1]
            # ys[sharp_turns] = ys[sharp_turns-1]
            xs[sharp_turns] = (xs[sharp_turns-1] + xs[sharp_turns+1]) / 2
            ys[sharp_turns] = (ys[sharp_turns-1] + ys[sharp_turns+1]) / 2
            n += 1

    if n == nmax:
        print(f'warning: smooth_bank did not converge in {n} steps\n')
    else:
        print(f'smooth_bank converged in {n} steps\n')
        
    return xs, ys

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

    return S

if __name__ == "__main__":
    # tif_fname = r'/sciclone/data10/wangzg/DEM/npz/sc_ll_7.npz'  # directory of DEM data
    tif_fname = 'LA_CONED_merged_3_UTM15_5m.tif'
    search_steps = 100
    river_threshold = np.array([35, 500])
    
    # %%
    if pathlib.Path(tif_fname).suffix == ".npz":
        S = loadz(tif_fname)
    elif pathlib.Path(tif_fname).suffix == ".tif" : 
        S = Tif2XYZ(tif_fname=tif_fname)
    else:
        raise Exception("Unknown DEM format.")

    # %%
    xyz, l2g, curv = get_all_points_from_shp('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_test.shp')

    range_arcs = []
    bank_arcs = []
    improved_thalweg_arcs = []
    search_length = river_threshold[1] * 1.1

    for i, idx in enumerate(l2g):
        print(f'Arc {i+1} of {len(l2g)}')
        
        line = xyz[idx, :]

        line_cplx = np.squeeze(line.view(np.complex128))
        angles = np.angle(np.diff(line_cplx))
        angle_diff0 = np.diff(angles)
        angle_diff = np.diff(angles)
        angle_diff[angle_diff0 > np.pi] -= 2 * np.pi
        angle_diff[angle_diff0 < -np.pi] += 2 * np.pi
        perp = angles[:-1] + angle_diff / 2 - np.pi / 2
        perp = np.r_[angles[0] - np.pi / 2, perp, angles[-1] - np.pi / 2]

        # %%
        xt_right = line[:, 0] + search_length * np.cos(perp)
        yt_right = line[:, 1] + search_length * np.sin(perp)

        xt_left = line[:, 0] + search_length * np.cos(perp + np.pi)
        yt_left = line[:, 1] + search_length * np.sin(perp + np.pi)

        range_arcs += [SMS_ARC(points=np.c_[xt_left, yt_left]), SMS_ARC(points=np.c_[xt_right, yt_right])]

        # %%
        dx = S.lon[1] - S.lon[0]
        dy = S.lat[1] - S.lat[0]

        # %%
        x = line[:, 0]
        y = line[:, 1]
        # %%
        x, y = improve_thalwegs(S, x, y, xt_left, yt_left, xt_right, yt_right, search_steps=int(0.2*search_steps))
        improved_thalweg_arcs += [SMS_ARC(points=np.c_[x, y])]

        thalweg_eta = np.maximum(0.0, (y - 3313760.0))/(3367300.0 - 3313760.0) * 8

        x_banks_left, y_banks_left = get_bank(S, x, y, thalweg_eta, xt_left, yt_left, search_steps=search_steps)
        x_banks_left, y_banks_left = nudge_bank(line, perp+np.pi, x_banks_left, y_banks_left, dist=river_threshold)
        x_banks_left, y_banks_left = smooth_bank(x_banks_left, y_banks_left)

        x_banks_right, y_banks_right = get_bank(S, x, y, thalweg_eta, xt_right, yt_right, search_steps=search_steps)
        x_banks_right, y_banks_right = nudge_bank(line, perp, x_banks_right, y_banks_right, dist=river_threshold)
        x_banks_right, y_banks_right = smooth_bank(x_banks_right, y_banks_right)


        bank_arcs.append(SMS_ARC(points=np.c_[x_banks_left, y_banks_left]))
        bank_arcs.append(SMS_ARC(points=np.c_[x_banks_right, y_banks_right]))

    SMS_MAP(arcs=improved_thalweg_arcs).writer(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_improved_thalwegs.map')
    map = SMS_MAP(arcs=range_arcs)
    map.writer(filename=f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_test_range.map')

    map = SMS_MAP(arcs=bank_arcs)
    map.writer(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/LA_test_bank.map')

    # %%
    pass
