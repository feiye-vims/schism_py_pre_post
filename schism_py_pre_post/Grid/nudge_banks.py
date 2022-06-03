import numpy as np
import shapefile
from scipy import spatial
from schism_py_pre_post.Grid.SMS import get_all_points_from_shp, replace_shp_pts, find_pts_in_shpfiles


if __name__ == '__main__':
    state = 'NH'
    # original
    river_bank_shpfname = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/{state}_bank.shp'
    # based on redistributed arc pts
    thalweg_buffer_poly_shpfname = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/2_{state}_stream_35_buffer_poly_redist_50m.shp'
    # outputs
    river_bank_outfname = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/3_{state}_river_bank_nudged.shp'

    river_banks_pts, shape_pts_l2g, _ = get_all_points_from_shp(river_bank_shpfname)
    thalweg_buffer_pts, *dummy = get_all_points_from_shp(thalweg_buffer_poly_shpfname)

    mindist_bank_idx = spatial.cKDTree(thalweg_buffer_pts).query(river_banks_pts)[1]

    too_close_idx = find_pts_in_shpfiles([thalweg_buffer_poly_shpfname], pts=river_banks_pts)
    river_banks_pts[too_close_idx] = thalweg_buffer_pts[mindist_bank_idx[too_close_idx]]

    replace_shp_pts(
        river_bank_shpfname,  # original
        river_banks_pts, shape_pts_l2g,  # updated pts info
        river_bank_outfname,
    )

    # Optional
    # np.savetxt(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/{state}_river_bank_redist_nudged.xyz', river_banks_pts)

    # w = shapefile.Writer(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/{state}/{state}_river_bank_redist_nudged.shp')
    # w.field('name', 'C')
    # w.line(river_banks_pts[l2g])
    # w.record('linestring1')
    # w.close()
