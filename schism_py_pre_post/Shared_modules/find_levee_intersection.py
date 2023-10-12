import geopandas as gpd
from itertools import combinations

# read levee shapefile
levee_shp_fname = '/sciclone/schism10/Hgrid_projects/STOFS3D-V6/v18-s2_v1/levee_v3.shp'
gdf1 = gpd.read_file(levee_shp_fname)

# read other shapefile
other_shp_fname = '/sciclone/schism10/Hgrid_projects/STOFS3D-V6/v18-s2_v1/coastal_v43-s2_v1.shp'
gdf2 = gpd.read_file(other_shp_fname)

# find intersection points 
intersections = []

# Iterate over the other GeoDataFrame
sindex = gdf2.sindex
intersection_points = []
for i, line in enumerate(gdf1['geometry']):
    # Find approximate matches with the bounding box (this is fast)
    possible_matches_index = list(sindex.intersection(line.bounds))
    if not possible_matches_index:
        continue

    possible_matches = gdf2.iloc[possible_matches_index]
    
    # Find actual intersections
    # precise_matches = possible_matches[possible_matches.intersects(line)]
    
    for match_line in possible_matches['geometry']:
        intersection = line.intersection(match_line)
        
        if intersection.geom_type == 'Point':
            intersection_points.append(intersection)
        elif intersection.geom_type == 'MultiPoint':
            intersection_points.extend(list(intersection.geoms))

# Convert to a GeoDataFrame
intersection_gdf = gpd.GeoDataFrame({'geometry': intersection_points})
# write to shapefile
intersection_gdf.to_file('/sciclone/schism10/Hgrid_projects/STOFS3D-V6/v18-s2_v1/levee_v3_intersection.shp')

print(intersection_gdf)