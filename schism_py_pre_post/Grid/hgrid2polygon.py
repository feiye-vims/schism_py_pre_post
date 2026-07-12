import geopandas as gpd
from shapely.ops import polygonize, unary_union
from pylib import read


hg = read("/sciclone/schism10/feiye/STOFS3D-v8/R29i1/hgrid.gr3")
hg.write_shp("/sciclone/schism10/feiye/STOFS3D-v8/R29i1/hgrid.shp", fmt=2)

# convert to polygon
shp = gpd.read_file("/sciclone/schism10/feiye/STOFS3D-v8/R29i1/hgrid.shp")
lines = unary_union(shp.geometry)
polys = list(polygonize(lines))

domain_poly = unary_union(polys)

out = gpd.GeoDataFrame(
    {"name": ["hgrid_domain"]},
    geometry=[domain_poly],
    crs=shp.crs
)

out.to_file(
    "/sciclone/schism10/feiye/STOFS3D-v8/R29i1/hgrid_polygon.gpkg",
    driver="GPKG"
)
