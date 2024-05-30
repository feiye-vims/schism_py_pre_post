import geopandas as gpd

nwm_hydrofabric = gpd.read_file('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2v1_RiverMapper/shapefiles/LA_nwm_v1p2.shp')

# remove order 1 streams
nwm_hydrofabric = nwm_hydrofabric[nwm_hydrofabric['order_'] > 1]

nwm_hydrofabric.to_file('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2v1_RiverMapper/shapefiles/LA_nwm_v1p2_order2.gpkg', driver='GPKG')

pass