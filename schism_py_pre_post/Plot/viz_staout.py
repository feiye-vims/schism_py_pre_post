from pylib_essentials.schism_file import TimeHistory
from schism_py_pre_post.Grid.Bpfile import Bpfile
from matplotlib import pyplot as plt
from pathlib import Path

runs = [
    # '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R03a',
    # '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R03',
    # '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R02',
    # '/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/RUN24/',
    '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R12',
    '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R14',
]
line_styles = ['-', '--', '.-', '-.', ':', '-.', '--']
line_colors = ['r', 'b', 'g', 'k', 'c', 'm', 'y', '#FFA500', '#800080', '#00FFFF', '#008000', '#FF00FF', '#808080', '#800000', '#000080', '#FFFF00', '#00FF00', '#FF0000', '#0000FF', '#C0C0C0', '#008080', '#000000']

fig, ax = plt.subplots()

for i, [run, line_style] in enumerate(zip(runs, line_styles)):
    fname = f'{run}/outputs/staout_1'
    station_in = f'{run}/station.in'
    stations = Bpfile(station_in, cols=5)
    ts = TimeHistory.from_file(fname, columns=stations.st_id, start_time_str="2021-08-01 00:00:00")
    ts.df.iloc[:229,:160:32].plot(color=line_colors, style=line_style, ax=ax, label=Path(run).stem)

def onpick2(event):
    print('onpick2 line:', event.pickx, event.picky)
fig.canvas.mpl_connect('pick_event', onpick2)

plt.show()
