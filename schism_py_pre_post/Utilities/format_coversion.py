"""
Convert commonly used formats to each other
"""


from pathlib import Path
import numpy as np
import geopandas as gpd
from pylib import schism_bpfile
from pylib import read, save, read_schism_reg
from shapely.geometry import Polygon


def shp2reg(shp_fname: Path, output_epsg: int = 4326):
    """
    Convert a shapefile to a SCHISM region file

    Parameters:
    - shp_fname: Path, path to the shapefile, which has one or more polygons

    Outputs:
    - Writes SCHISM region files to the same directory as the shapefile,
      each named with the polygon index.
    """
    gdf = gpd.read_file(shp_fname).to_crs(epsg=output_epsg)

    for i, row in gdf.iterrows():
        # test if the row is a polygon
        if row['geometry'].geom_type == 'Polygon':
            bp = schism_bpfile(x=row['geometry'].exterior.xy[0], y=row['geometry'].exterior.xy[1])
            bp.write(f'{shp_fname.parent}/{shp_fname.stem}_{i}.reg')
        else:
            print(f'Row {i} is not a polygon')


def sample_usage_shp2reg():
    """
    Example usage of shp2reg function
    """            
    shp_fname = Path(
        '/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/'
        'STOFS-3D-Atl-shadow-VIMS/Pre_processing/Gr3/Drag/Shapefiles/'
        'drag_reduce_Oyster_Landing.shp'
    )
    shp2reg(shp_fname, output_epsg=4326)


def reg2shp(reg_fname: Path, output_epsg: int = 4326):
    """
    Convert a SCHISM region file to a shapefile

    Parameters:
    - reg_fname: Path, path to the SCHISM region file

    Outputs:
    - Writes a shapefile to the same directory as the region file
    """
    region = read_schism_reg(reg_fname)
    polygon = Polygon(np.r_[np.c_[region.x, region.y], np.c_[region.x[0], region.y[0]]])  # close the polygon
    gdf = gpd.GeoDataFrame({'geometry': [polygon]})
    gdf.to_file(reg_fname.with_suffix('.shp'))


def sample_usage_reg2shp():
    """
    Example usage of reg2shp function
    """
    reg_fname = Path(
        '/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/'
        'STOFS-3D-Atl-shadow-VIMS/Pre_processing/VDatum/chea_del_bay.reg'
    )
    reg2shp(reg_fname, output_epsg=4326)


if __name__ == "__main__":
    sample_usage_reg2shp()
    print('Done.')