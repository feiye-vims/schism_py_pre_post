from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Bpfile import Bpfile
import matplotlib.pyplot as plt


# nodes = [
#   [-(91+11/60+29.6/3600), 30+26/60+44.4/3600],
#   [-(90+8/60+10/3600), 29+56/60+5/3600]
# ]
# mybp = Bpfile(nodes=nodes)
# mybp.writer(out_file_name='/sciclone/schism10/feiye/ICOGS/BPfiles/missi_gages.bp')

conv_factor = 3.28084

myth = TimeHistory(file_name='/sciclone/schism10/feiye/ICOGS/RUN20a/PostP/elev.dat.missi_gages.RUN20a', sec_per_time_unit=86400, start_time_str='2005-08-01 00:00:00')
plt.plot(myth.df['datetime'], myth.df.iloc[:, 1:]*conv_factor)
plt.show()

pass
