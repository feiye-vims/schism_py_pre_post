from pylib_essentials.schism_file import TimeHistory
from schism_py_pre_post.Grid.Bpfile import Bpfile
from matplotlib import pyplot as plt
from pathlib import Path

runs = [
    '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R00a',
    '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R00b',
]
line_styles = ['-', '--']
line_colors = ['r', 'b', 'g', 'k', 'c', 'm', 'y']

for i, [run, line_style] in enumerate(zip(runs, line_styles)):
    fname = f'{run}/outputs/staout_1'
    station_in = f'{run}/station.in'
    stations = Bpfile(station_in, cols=5)
    ts = TimeHistory.from_file(fname, columns=stations.st_id, start_time_str="2018-08-01 00:00:00")
    if i == 0:
        ax = ts.df.iloc[:,::4].plot(color=line_colors, style=line_style, label=Path(run).stem)
    else:
        ts.df.iloc[:,::4].plot(color=line_colors, style=line_style, ax=ax, label=Path(run).stem)

pass
plt.show()
