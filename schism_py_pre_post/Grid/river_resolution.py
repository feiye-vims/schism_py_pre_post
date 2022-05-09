import numpy as np
import shapefile
from scipy import spatial


def get_all_points_from_shp(fname):
    sf = shapefile.Reader(fname)
    shapes = sf.shapes()
    xyz = np.empty((0, 2), dtype=float)
    for i, shp in enumerate(shapes):
        print(f'shp {i} of {len(shapes)}')
        xyz = np.append(xyz, shp.points, axis=0)

    return xyz
            

if __name__ == '__main__':
    river_banks_pts = get_all_points_from_shp('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/ME/ME_river_bank.shp')
    thalweg_pts = get_all_points_from_shp('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/ME/ME_river_stream.shp')

    # nearest neighbor
    tree = spatial.cKDTree(river_banks_pts)
    river_widths = tree.query(thalweg_pts)[0]

    # set limits on resolution
    river_resolution = np.maximum(50, np.minimum(500, river_widths*5))
    np.savetxt('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/ME/ME_river_resolution.xyz', np.c_[thalweg_pts, river_resolution])
