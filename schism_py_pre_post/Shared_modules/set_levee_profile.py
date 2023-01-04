import shapefile
import numpy as np
from pylib import schism_grid, inside_polygon, proj
from schism_py_pre_post.Shared_modules.hotstart_proc import nearest_neighbour
from schism_py_pre_post.Download.download_nld import nld2map
import os


def set_levee_profile(gd=None, wdir='./'):

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()

    levee_names = ['LA_levees', 'FL_levees']
    levee_name_str = "_".join(levee_names)
    levee_xyz = np.zeros((0, 3), dtype=float)

    for levee_name in levee_names:
        # read levee heights as xyz
        _, xyz = nld2map(nld_fname=f'{levee_info_dir}/{levee_name}/System.geojson')
        levee_xyz = np.r_[levee_xyz, xyz]
    levee_x = levee_xyz[:, 0]
    levee_y = levee_xyz[:, 1]
    levee_height = levee_xyz[:, 2]
    levee_height[levee_height < 1] = 27  # raise very low levees to 9 meters
    levee_height *= 0.3048  # convert to meters
    # plt.plot(np.sort(levee_height))
    # plt.show()

    if gd is None:
        gd = schism_grid(f'{wdir}/hgrid.ll')  # ; gd.save(f'{wdir}/hgrid.pkl')

    gd.lon = gd.x
    gd.lat = gd.y
    gd.proj(prj0='epsg:4326', prj1='esri:102008')  # this overwrites gd.x, gd.y

    # find levee center line points in hgrid, use UTM to avoid truncation error
    shapefile_names = [
        f"{levee_info_dir}/Polygons/la_levee_center_line_buffer_102008.shp",
        f"{levee_info_dir}/Polygons/fl_levees_buffer_10m_102008.shp",
    ]
    ilevee = np.zeros(gd.dp.shape)
    for shapefile_name in shapefile_names:
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()
        for i, shp in enumerate(shapes):
            print(f'shp {i+1} of {len(shapes)}')
            poly_xy = np.array(shp.points).T
            ilevee += inside_polygon(np.c_[gd.x, gd.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
    ilevee = ilevee.astype('bool')

    gd.save(f'{wdir}/{levee_name_str}.gr3', value=ilevee)

    II = nearest_neighbour(np.c_[gd.lon[ilevee], gd.lat[ilevee]], np.c_[levee_x, levee_y])
    dist = np.sqrt((gd.lon[ilevee] - levee_x[II])**2 + (gd.lat[ilevee] - levee_y[II])**2)
    short_dist = dist < 0.01  # degree, roughly 1000 m

    # gd.dp[:] = 0
    idx_levee_in_range = np.argwhere(ilevee)[:, 0][short_dist]
    gd.dp[idx_levee_in_range] = - levee_height.astype(float)[II][short_dist]

    gd.x = gd.lon
    gd.y = gd.lat

    # os.system(f"cp {wdir}/hgrid_{levee_name_str}_loaded_ll.gr3 {wdir}/hgrid.ll")
    # proj(
    #     f'{wdir}/hgrid.ll', 0, 'epsg:4326',
    #     f'{wdir}/hgrid.utm.gr3', 0, 'epsg:26918',
    # )

    return gd  # levee loaded hgrid.ll

if __name__ == "__main__":
    '''
    Inputs under wdir, hgrid.ll
    Outputs to wdir: levee-loaded hgrid.ll
    '''

    wdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ/Test_levee_heights/'
    gd = set_levee_profile(gd=None, wdir=wdir)
    gd.write_hgrid(f'{wdir}/hgrid_levee_profile_loaded_ll.gr3')

    pass
