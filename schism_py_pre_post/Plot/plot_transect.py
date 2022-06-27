from cProfile import label
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
import numpy as np
import matplotlib.pyplot as plt


runids = ['RUN20', 'RUN21f']

start_strs = ['2021-08-24 00:00:00', '2021-05-01 00:00:00']
sss = []
for runid, start_str in zip(runids, start_strs):
    ss = TimeHistory(file_name=f'/sciclone/schism10/feiye/ICOGS/{runid}/vsource.th',
                     start_time_str=start_str)
    mask = (ss.df['datetime'] > '2021-8-1') & (ss.df['datetime'] <= '2021-9-10')
    ss.df = ss.df.loc[mask]
    sss.append(ss.df.max(axis=1))

plt.plot(ss.df['datetime'], ss.df[])

np.c_[sss[0], sss[1]]

runids = ['RUN10a', 'RUN10b', 'RUN10d', 'RUN10e', 'RUN10g']
markers = ['o', '*', '+' , 's', 'd']
ths = []
for runid in runids:
    th = TimeHistory(file_name=f'/sciclone/schism10/feiye/ICOGS/{runid}/PostP/elev.dat.missi.{runid}',
                     start_time_str='2021-05-01 00:00:00', sec_per_time_unit=86400)
    mask = (th.df['datetime'] > '2021-9-2') & (th.df['datetime'] <= '2021-9-4')
    th.df = th.df.loc[mask]
    th.df = th.df.append(th.df.max(axis=0), ignore_index=True)
    ths.append(th)

for th, runid, marker in zip(ths, runids, markers):
    th.df.iloc[-1, 1:].to_csv(f'/sciclone/schism10/feiye/ICOGS/{runid}/PostP/elev.dat.missi.{runid}.tavg.txt')
    plt.plot(th.df.iloc[-1, 1:], label=runid, marker=marker)
    # plt.plot(th.df.iloc[0, 1:], label=runid)
plt.legend()
plt.show()
# np.transpose(np.r_[[ths[0].df.iloc[-1, 1:].values], [ths[1].df.iloc[-1, 1:].values]])

pass