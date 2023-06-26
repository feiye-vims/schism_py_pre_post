from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Bpfile import Bpfile
from matplotlib import pyplot as plt


fname = '/sciclone/schism10/feiye/STOFS3D-v6/Runs/RUN15/outputs/staout_1'
station_in = '/sciclone/schism10/feiye/STOFS3D-v6/Runs/RUN15/station.in'

stations = Bpfile(station_in, cols=5)
ts = TimeHistory(fname, columns=stations.st_id, start_time_str="2017-12-01 00:00:00")
ts.df.set_index('datetime').iloc[:, 0:180:20].plot()
# plt.ylim([0, 38])
plt.show()

pass
