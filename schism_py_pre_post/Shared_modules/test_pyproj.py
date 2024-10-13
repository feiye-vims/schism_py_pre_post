"""Test pyproj module"""

import pyproj
from pyproj import Transformer

# Define the transformation pipeline string
nad832011_navd88geoid18_to_itrf14_2010_xgeoid20b = """+proj=pipeline
  +step +proj=axisswap +order=2,1
  +step +proj=unitconvert +xy_in=deg +xy_out=rad
  +step +proj=vgridshift +grids=us_noaa_g2018u0.tif +multiplier=1
  +step +proj=cart +ellps=GRS80
  +step +inv +proj=helmert +x=1.0053 +y=-1.9092 +z=-0.5416 +rx=0.0267814
        +ry=-0.0004203 +rz=0.0109321 +s=0.00037 +dx=0.0008 +dy=-0.0006
        +dz=-0.0014 +drx=6.67e-05 +dry=-0.0007574 +drz=-5.13e-05 +ds=-7e-05
        +t_epoch=2010 +convention=coordinate_frame
  +step +inv +proj=cart +ellps=GRS80
  +step +proj=vgridshift +grids=/sciclone/schism10/Hgrid_projects/DEMs/xGEOID20B.tif +multiplier=-1
  +step +proj=unitconvert +xy_in=rad +xy_out=deg
  +step +proj=axisswap +order=2,1"""

# Create the transformer using pyproj
transformer = Transformer.from_pipeline(nad832011_navd88geoid18_to_itrf14_2010_xgeoid20b)

# lon, lat, z = 35.685591, -75.447685, 19.00885703
# x, y, z_trans = transformer.transform(lon, lat, z)

x, y, z_trans, _ = transformer.transform(38.88950026013818, -77.03529973563386, 0, 2010.0)

# Output the transformed coordinates
print(f"Original coordinates (lon, lat): ({lon}, {lat}, {z})")
print(f"Transformed coordinates (xGEOID20B): ({x}, {y}, {z_trans})")
