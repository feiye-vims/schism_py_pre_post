import geopandas as gpd
from shapely.ops import unary_union, snap
from joblib import Parallel, delayed

wdir = '/sciclone/schism10/feiye/Test/'

# Load the shapefile
gdf = gpd.read_file(f'{wdir}/test.shp')

# Set a small snapping tolerance
tolerance = 1e-5  # about 1 meter in degrees

# Function to process each polygon
def process_polygon(polygon, others):
    # Snap the polygon to the unary union of all other polygons
    snapped_polygon = snap(polygon, unary_union(others), tolerance)
    return snapped_polygon

# Parallel processing function
def process_in_parallel(gdf):
    processed_polygons = Parallel(n_jobs=-1)(delayed(process_polygon)(row.geometry, gdf.geometry.drop(index))
                                             for index, row in gdf.iterrows())
    return processed_polygons

# Run the process
gdf['geometry'] = process_in_parallel(gdf)

# Save the corrected shapefile
gdf.to_file('fixed_shapefile.shp')
