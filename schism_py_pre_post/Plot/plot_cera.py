# import the Cartopy/Matplotlib libraries
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import cartopy.crs as ccrs

# libraries for data handling
import numpy as np
import netCDF4

#For the purposes of this tutorial, we'll be ignoring warnings from functions
import warnings
warnings.filterwarnings('ignore')

wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/CERA_salinity_map/'

# Open the NetCDF file
# https://noaa-nos-stofs3d-pds.s3.amazonaws.com/index.html#STOFS-3D-Atl/para/stofs_3d_atl.20240506/
# s3://noaa-nos-stofs3d-pds/STOFS-3D-Atl/para/stofs_3d_atl.20240506

# reading in the NetCDF file
myfile = netCDF4.Dataset(f'{wdir}/schout_adcirc_20240506.nc')
#print(myfile)

# getting all variables names as contained in the file by printing the dictionary ’keys’
all_variable_names = myfile.variables.keys()
print(all_variable_names)

# getting a list of all variable attributes (metadata) in the file
# note, this does not print any data but just the attributes
#for attrs in myfile.variables.values():
  # print(attrs)

# accessing the attributes (metadata) of the variable ’x’ (longitudes) and  the variable ’y’ (latitudes)
print(myfile.variables['x'], myfile.variables['y'], myfile.variables['element'])

# Open the NetCDF file
# https://noaa-nos-stofs3d-pds.s3.amazonaws.com/index.html#STOFS-3D-Atl/para/stofs_3d_atl.20240506/
# s3://noaa-nos-stofs3d-pds/STOFS-3D-Atl/para/stofs_3d_atl.20240506

# reading in the NetCDF file
myfile = netCDF4.Dataset(f'{wdir}/stofs_3d_atl.t12z.field2d_f001_012.nc')
#print(myfile)

# getting all variables names as contained in the file by printing the dictionary ’keys’
all_variable_names = myfile.variables.keys()
print(all_variable_names)

# getting a list of all variable attributes (metadata) in the file
# note, this does not print any data but just the attributes
#for attrs in myfile.variables.values():
  # print(attrs)

# accessing the attributes (metadata) of the variable ’x’ (longitudes) and  the variable ’y’ (latitudes)
print(myfile.variables['SCHISM_hgrid_node_x'], myfile.variables['SCHISM_hgrid_node_y'], myfile.variables['SCHISM_hgrid_face_nodes'])


# The plt.figure object is the Matplotlib container that holds all inner elements, including the plotting area (map).
plt.figure(figsize=(14,7)) # width, height (inches)

# set the map projection and map extent
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-86, -97, 26, 33])

# Open the NetCDF file
# https://noaa-nos-stofs3d-pds.s3.amazonaws.com/index.html#STOFS-3D-Atl/para/stofs_3d_atl.20240506/
# s3://noaa-nos-stofs3d-pds/STOFS-3D-Atl/para/stofs_3d_atl.20240506

myfile = f'{wdir}/schout_adcirc_20240506.nc'
# myfile = '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA2/extract/schout_adcirc_1.nc'
vars = netCDF4.Dataset(myfile).variables

# access grid variables (x, y, element)
x = vars['x'][:]
y = vars['y'][:]
elems = vars['element'][:, :] - 1  # Move to 0-indexing by subtracting 1, elements indexing starts with '1' in netcdf file

# read the 'attribute name' data array from the input file
data = vars['zeta_max'][:]

 # matplotlib: triangulation
triang = tri.Triangulation(x, y, triangles=elems)

# check if data array is masked
if data.mask.any():
  # -99999 entries in 'data' array are usually masked, mask all corresponding triangles
  point_mask_indices = np.where(data.mask)
  tri_mask = np.any(np.in1d(elems, point_mask_indices).reshape(-1, 3), axis=1)
  triang.set_mask(tri_mask)

levels = np.linspace(0, 5, num=30)

# Matplotlib: Color plot
c = plt.tricontourf(triang, data, levels=levels, cmap=plt.cm.jet, extend='both')
c.cmap.set_under("#000066")
c.cmap.set_over("#880066")
plt.colorbar(c, ticks=levels, pad=0.02)

# add a map title
ax.set_title('Maximum Water Levels [m]')

plt.show()

#plt.savefig('cera_max_water_levels.png',bbox_inches = "tight", format='png', dpi=300) 


# The plt.figure object is the Matplotlib container that holds all inner elements, including the plotting area (map).
plt.figure(figsize=(14,7)) # width, height (inches)

# set the map projection and map extent
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-86, -97, 26, 33])

# Open the NetCDF file
# https://noaa-nos-stofs3d-pds.s3.amazonaws.com/index.html#STOFS-3D-Atl/para/stofs_3d_atl.20240506/
# s3://noaa-nos-stofs3d-pds/STOFS-3D-Atl/para/stofs_3d_atl.20240506

myfile = f"{wdir}/stofs_3d_atl.t12z.field2d_f001_012.nc"
myfile = "/sciclone/home/feiye/TEMP/v2.1/stofs_3d_atl.t12z.field2d_f025_036.nc"
# myfile = '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA2/extract/schout_UV4.5m_1.nc'
vars = netCDF4.Dataset(myfile).variables

# access grid variables (x, y, element)
x = vars['SCHISM_hgrid_node_x'][:]
y = vars['SCHISM_hgrid_node_y'][:]
elems = vars['SCHISM_hgrid_face_nodes'][:, :] - 1  # Move to 0-indexing by subtracting 1, elements indexing starts with '1' in netcdf file

# read the 'attribute name' data array from the input file
mask = vars['temp_surface'][:].mask
data = vars['temp_surface'][:]

 # matplotlib: triangulation
triang = tri.Triangulation(x, y, triangles=elems)

# check if data array is masked
if mask.any():
  # -99999 entries in 'data' array are usually masked, mask all corresponding triangles
  point_mask_indices = np.where(mask)
  tri_mask = np.any(np.in1d(elems, point_mask_indices).reshape(-1, 3), axis=1)
  triang.set_mask(tri_mask)

levels = np.linspace(0, 30, num=30)

# Matplotlib: Color plot
c = plt.tricontourf(triang, data[0, :], levels=levels, cmap=plt.cm.jet, extend='both')
c.cmap.set_under("#000066")
c.cmap.set_over("#880066")
plt.colorbar(c, ticks=levels, pad=0.02)

# add a map title
ax.set_title('Surface Salinity [psu]')

plt.show()

#plt.savefig('cera_max_water_levels.png',bbox_inches = "tight", format='png', dpi=300)
pass