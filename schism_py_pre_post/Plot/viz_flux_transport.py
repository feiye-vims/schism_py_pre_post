# %%
import numpy as np
import matplotlib.pyplot as plt
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory


flux = TimeHistory('/sciclone/scr10/feiye/STOFS3D-v4/RUN23p9/outputs/flux.out',
                   sec_per_time_unit=86400)

# %%
plt.plot(flux.df['datetime'], flux.df[1])
plt.xticks(rotation=20)
pass
# %%
