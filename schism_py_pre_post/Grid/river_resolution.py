import numpy as np
import shapefile
from scipy import spatial
from schism_py_pre_post.Grid.SMS import get_all_points_from_shp


if __name__ == '__main__':
    state = 'NH'
    # river_bank_fname = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/3_{state}_river_bank_nudged.shp'
    river_bank_fname = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/{state}_bank.shp'
    thalweg_fname = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/{state}_stream_redist_50m.shp'

    outfile = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/{state}_river_resolution.xyz'

    thalweg_pts, l2g, curv = get_all_points_from_shp(thalweg_fname)
    river_banks_pts, *dummy = get_all_points_from_shp(river_bank_fname)

    # nearest neighbor
    tree = spatial.cKDTree(river_banks_pts)
    river_widths = tree.query(thalweg_pts)[0]

    # set limits on resolution
    # river_resolution = np.maximum(60, np.minimum(400, river_widths*5))
    # river_resolution = np.maximum(30, np.minimum(400, river_widths*4))  # MSAL 
    river_resolution = np.maximum(60, np.minimum(400, river_widths*10))  # ME

    R = 1.0/(curv+1e-10)
    W = np.maximum(30, np.minimum(400, river_widths*5))
    river_resolution = np.minimum(0.5 * R, W)
    river_resolution = np.maximum(30, np.minimum(400, river_resolution))

    np.savetxt(outfile, np.c_[thalweg_pts, river_resolution])
    # np.savetxt(outfile1, np.c_[thalweg_pts, curv])
