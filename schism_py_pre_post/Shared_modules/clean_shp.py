'''
Clean up geometry of shapefile
'''


import geopandas as gpd


gdf = gpd.read_file('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v37s2_RiverMapper/shapefiles/rivers_v37.3.shp')

# Add original index as a new column
gdf['orig_index'] = gdf.index

# Perform intersection to find overlaps
overlap_gdf = gpd.overlay(gdf, gdf, how='intersection', keep_geom_type=False)

# Remove self-intersections (where a polygon intersects itself)
overlap_gdf = overlap_gdf[overlap_gdf['orig_index_1'] != overlap_gdf['orig_index_2']]

# Add area column
overlap_gdf['area'] = overlap_gdf.area

# Sort by index pairs and area, and retain only the largest polygon for each overlap
overlap_gdf = overlap_gdf.sort_values(by=['orig_index_1', 'orig_index_2', 'area'], ascending=[True, True, False])
largest_overlap_gdf = overlap_gdf.drop_duplicates(subset=['orig_index_1', 'orig_index_2'], keep='first')

largest_overlap_gdf = largest_overlap_gdf.explode(index_parts=True).reset_index(drop=True)

# Visualize or save the largest polygons
largest_overlap_gdf.plot()  # Visualize