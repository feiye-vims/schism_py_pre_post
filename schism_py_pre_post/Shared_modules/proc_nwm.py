import geopandas as gpd

nwm_hydrofabric = gpd.read_file('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v46/Shapefiles/ecgc.shp')

# remove order 1 streams
nwm_hydrofabric = nwm_hydrofabric[nwm_hydrofabric['order_'] > 1]

nwm_hydrofabric.to_file('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v46/Shapefiles/ecgc_order2.shp')

pass