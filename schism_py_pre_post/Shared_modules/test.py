from pylib import harmonic_analysis
import numpy as np

elev = np.loadtxt('/sciclone/data10/feiye/vims20/work/ChesBay/RUN200p/elev.dat.more.york2_moved')
amps = []
for i in range(1, len(elev[0])):
    ha = harmonic_analysis(elev[:, i], 1/24, tidal_names=['M2'])
    amps.append(ha.amplitude[1])

# %%
import pygrib
import copy
import numpy as np

filename = '/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/estofs.t00z.conus.east.cwl.grib2'
outfile = '/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/test.grib2'

gr = pygrib.open(filename)
with open(outfile, 'wb') as outgrb:
    for g in gr:
        # print(g.typeOfLevel, g.level, g.name, g.validDate, g.analDate, g.forecastTime)
        # vals = copy.deepcopy(g['values'])
        # vals[vals.mask == False] = 0
        g['values'] -= g['values']
        g['values'] -= g['values']
        g['values'] = g['values'] + 10.0
        msg = g.tostring()
        outgrb.write(msg)
outgrb.close()
gr.close()
pass

# %%
import xarray as xr
from cfgrib.xarray_to_grib import to_grib


filename = '/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/estofs.t00z.conus.east.cwl.grib2'
ds=xr.open_dataset(filename, engine='cfgrib')

ds.to_netcdf('test.nc')
 

# %%
from pyproj import Proj, transform
inProj = Proj('esri:102009')
outProj = Proj('epsg:4326')

x2,y2 = transform(inProj,outProj, 2681.9,-263.8)

pass
