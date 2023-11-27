from pylib_essentials.schism_file import TimeHistory
from schism_py_pre_post.Grid.Bpfile import Bpfile
from matplotlib import pyplot as plt
from pathlib import Path

runs = [
    '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R00a',
    '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R01',
    #'/sciclone/schism10/feiye/STOFS3D-v7/Runs/R02',
]
line_styles = ['-', '--', '.-', '-.', ':', '-.', '--']
line_colors = ['r', 'b', 'g', 'k', 'c', 'm', 'y']

fig, ax = plt.subplots()

for i, [run, line_style] in enumerate(zip(runs, line_styles)):
    fname = f'{run}/outputs/staout_1'
    station_in = f'{run}/station.in'
    stations = Bpfile(station_in, cols=5)
    ts = TimeHistory.from_file(fname, columns=stations.st_id, start_time_str="2018-08-01 00:00:00")
    ts.df.iloc[:,:].plot(color=line_colors, style=line_style, ax=ax, label=Path(run).stem)

def onpick2(event):
    print('onpick2 line:', event.pickx, event.picky)
fig.canvas.mpl_connect('pick_event', onpick2)

plt.show()
