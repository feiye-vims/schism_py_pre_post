"Visualize NWM time series at featureIDs"
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

ds = xr.open_mfdataset('/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/CHRTOUT/2019/201906*CHRTOUT*')

pass

featureIDs = [19269176]  # Mississippi
varname = 'streamflow'
for featureID in featureIDs:
    data = ds[varname][ds['feature_id'] == featureID].values.flatten()

    plt.plot(time, data, label=f'featureID {featureID}')
plt.xlabel('Time')
plt.ylabel(varname)
plt.legend()
plt.title(f'NWM {varname} Time Series')
plt.show()

