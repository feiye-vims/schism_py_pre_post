# %%
import numpy as np
import matplotlib.pyplot as plt
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from datetime import datetime

labels = [
    # 'RUN23p7, AVISO-0.55m except for North of 40.5N',
    # 'RUN23p10, AVISO-0.5m',
    # 'RUN23p1, Original HYCOM',
    'R09d'
]

fluxes = [
    # TimeHistory('/sciclone/scr10/feiye/STOFS3D-v4/RUN23p7/outputs/flux.out',
    #             start_time_str="2021-05-01 00:00:00" , sec_per_time_unit=86400),
    # TimeHistory('/sciclone/scr10/feiye/STOFS3D-v4/RUN23p10/outputs/flux.out',
    #             start_time_str="2021-05-01 00:00:00" , sec_per_time_unit=86400),
    # TimeHistory('/sciclone/scr10/lcui01/ICOGS3D/outputs_RUN23p1/flux.out',
    #             start_time_str="2021-05-01 00:00:00" , sec_per_time_unit=86400),
    TimeHistory('/sciclone/scr10/feiye/R09d/outputs/flux.out',
                start_time_str="2024-03-05 00:00:00" , sec_per_time_unit=86400),
]

# %%
plt.figure(figsize=(12, 7))
for label, flux in zip(labels, fluxes):
    # II = (datetime(2021, 7, 5) <= np.array(flux.datetime)) * (np.array(flux.datetime) < datetime(2021, 8, 5))
    # print(np.mean(flux.df[1][II]))
    # plt.plot(flux.df['datetime'][II], flux.df[1][II], label=label)

    for col in flux.df.columns[1:]:
        plt.plot(flux.df['datetime'], flux.df[col], label=label + col)

plt.xticks(rotation=20)
plt.legend()
plt.show()
pass
# %%
