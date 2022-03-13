from posixpath import dirname
import numpy as np
from schism_py_pre_post.Grid.SourceSinkIn import source_sink
import matplotlib.pyplot as plt


ss = source_sink(
  source_dir='/sciclone/schism10/feiye/ICOGS/RUN22/Pumps/',
  start_time_str='2021-08-10 00:00:00'
)

mask = (ss.vsource.df['datetime'] > '2021-8-28') & (ss.vsource.df['datetime'] <= '2021-9-14')
df = ss.vsource.df[mask]
df_mean = df.mean(axis=0)
miss_idx = np.argmax(df_mean.values)
df_miss = df.iloc[:, miss_idx+1]  # 0 is time, so +1

# amplification = 1.25e6 * 0.028316847 / max(df_miss)
# ss.vsource.df.iloc[:, miss_idx+1] = ss.vsource.df.iloc[:, miss_idx+1] * amplification

offset = 9293.589185399998
ss.vsource.df.iloc[:, miss_idx+1] = ss.vsource.df.iloc[:, miss_idx+1] + offset

ss.vsource.writer(out_file_name='/sciclone/schism10/feiye/ICOGS/RUN22/Amp_source_at_mississippi/vsource.th')

pass