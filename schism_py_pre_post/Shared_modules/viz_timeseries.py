from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
import matplotlib.pyplot as plt
import numpy as np


vsource = TimeHistory(file_name=f'/sciclone/schism10/feiye/ICOGS/RUN22/Pumps/vsource.th', start_time_str='2021-08-10 00:00:00')
mask = (vsource.df['datetime'] >= '2021-8-10') & (vsource.df['datetime'] <= '2021-9-14')
df = vsource.df[mask]
df_mean = df.mean(axis=0)
miss_idx = np.argmax(df_mean.values)
df_miss = df.iloc[:, miss_idx+1]  # 0 is time, so +1
print(f"avg miss flow = {df_miss.mean()}")

plt.plot(df['datetime'], df_miss*35.31466621266132)
plt.show()
pass
