"""
Convert commonly used formats to each other
"""


from pathlib import Path
import geopandas as gpd
from pylib import schism_bpfile
from pylib import read, save


# read polygon shapefile
shp_fname = Path(
    '/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/'
    'STOFS-3D-Atl-shadow-VIMS/Pre_processing/Gr3/Drag/Shapefiles/'
    'drag_reduce_Oyster_Landing.shp'
)
gdf = gpd.read_file(shp_fname).to_crs(epsg=4326)

for i, row in gdf.iterrows():
    # test if the row is a polygon
    if row['geometry'].geom_type == 'Polygon':
        bp = schism_bpfile(x=row['geometry'].exterior.xy[0], y=row['geometry'].exterior.xy[1])
        bp.write(f'{shp_fname.parent}/{shp_fname.stem}_{i}.reg')
    else:
        print(f'Row {i} is not a polygon')



print('Done')