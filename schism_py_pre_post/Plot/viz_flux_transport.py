# %%
import numpy as np
import matplotlib.pyplot as plt
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from datetime import datetime


fluxes = [
    TimeHistory('/sciclone/scr10/feiye/STOFS3D-v4/RUN23p8/outputs/flux.out',
                start_time_str="2021-05-01 00:00:00" , sec_per_time_unit=86400),
    TimeHistory('/sciclone/scr10/feiye/STOFS3D-v4/RUN23p9/outputs/flux.out',
                start_time_str="2021-05-01 00:00:00" , sec_per_time_unit=86400),
    TimeHistory('/sciclone/scr10/lcui01/ICOGS3D/outputs_RUN23p1/flux.out',
                start_time_str="2021-05-01 00:00:00" , sec_per_time_unit=86400)


# %%
for flux in fluxes:
    II = (datetime(2021, 5, 5) <= np.array(flux.datetime)) * (np.array(flux.datetime) < datetime(2021, 5, 28))
    print(np.mean(flux.df[1][II]))
    plt.plot(flux.df['datetime'][II], flux.df[1][II])
plt.xticks(rotation=20)
pass
# %%
