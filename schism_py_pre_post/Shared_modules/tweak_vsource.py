import numpy as np
from pylib_experimental.schism_file import source_sink

wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I17b/Source_sink/'

ss = source_sink.from_ncfile(f'{wdir}/original_source_sink/source.nc', start_time_str='2021-08-01 00:00:00')

mask = (ss.vsource.df.index > '2021-8-28') & (ss.vsource.df.index <= '2021-9-14')
df = ss.vsource.df[mask]
df_mean = df.mean(axis=0)
mississippi_idx = np.argmax(df_mean.values)  # Mississippi river is the largest source in the domain

offset = 9293.589185399998
ss.vsource.df.iloc[:, mississippi_idx] = ss.vsource.df.iloc[:, mississippi_idx] + offset

ss.writer(f'{wdir}/adjusted_source_sink/')

pass
