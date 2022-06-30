import xarray
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.cm import ScalarMappable
import numpy as np


# data = xarray.open_dataset('/sciclone/home10/feiye/Sync/puertori.grib2', engine='cfgrib')
data = xarray.open_dataset('/sciclone/home10/feiye/Sync/conus.grib2', engine='cfgrib')

fig, ax = plt.subplots()
vmin = -1; vmax = 1
levels = np.linspace(vmin, vmax, 11)
C = ax.pcolor(data['longitude'], data['latitude'], data['unknown'][-1, :, :], vmin=vmin, vmax=vmax, cmap=cm.jet)
fig.colorbar(
   ScalarMappable(norm=C.norm, cmap=C.cmap),
   ticks=levels
)
plt.show()


fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(data['longitude'], data['latitude'], data['unknown'][0, :, :], cmap=cm.jet)
plt.show()

data.to_netcdf('netcdf_file.nc')