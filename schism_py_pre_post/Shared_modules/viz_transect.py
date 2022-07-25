#%%
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from pylib import schism_grid
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, rc

from schism_py_pre_post.Grid.Bpfile import Bpfile
import matplotlib

import plotly.io as pio
from IPython.display import set_matplotlib_formats


set_matplotlib_formats('svg')

matplotlib.rcParams['animation.embed_limit'] = 2**128
pio.renderers.default = 'iframe_connected'

class testAnimation:
    def __init__(self, datasets, runids, bp_x, bp_dp):
        self.plot_x = bp_x
        self.bp_x = bp_x
        self.bp_dp = bp_dp / 5  # scaled
        self.runids = runids

        self.datasets = datasets
        self.nt = datasets[0].shape[0]

        # First set up the figure, the axis, and the plot element we want to animate
        self.fig, ax = plt.subplots()
        ax.set_facecolor('white')
        plt.close()
        ax.set_xlim(( min(bp_x), max(bp_x)))
        ax.set_ylim((-10, 12))
        self.fig.set_size_inches(15, 15)

        # set lines (elevation from different runs)
        markers = ['o', '*', '+' , 's', 'd']
        self.lines = []
        for i, [runid, marker] in enumerate(zip(runids, markers)):
            line, = ax.plot([], [], lw=1, marker=marker, markevery=5, label=runid)
            self.lines.append(line)

        # set bottom (river bed)
        self.bottom, = ax.plot([], [], lw=2, label='river bed')
        ax.legend()

    # initialization function: plot the background of each frame
    def init(self):
        for i, _ in enumerate(self.runids):
            self.lines[i].set_data([], [])
        self.bottom.set_data([], [])
        return (self.lines+[self.bottom])

    # animation function. This is called sequentially  
    def animate(self, timestep):
        x = self.plot_x
        for k, data in enumerate(self.datasets):
            y = data[timestep, :]
            self.lines[k].set_data(x, y)
        self.bottom.set_data(self.bp_x, -self.bp_dp)
        return (self.lines+[self.bottom])

    def draw(self):
        global anim 
        anim = animation.FuncAnimation(self.fig, self.animate, init_func=self.init,
                                       frames=self.nt, interval=100, blit=True)

def get_transect_shape_from_grid(gridfile=None, bpfile=None):
    my_grid = schism_grid(gridfile)  # my_grid.save()
    my_bp = Bpfile(bpfile)

    bp_dp = my_grid.interp(pxy=np.c_[my_bp.nodes[:, 0], my_bp.nodes[:, 1]])
    bp_x, _ = np.unique(range(my_bp.n_nodes), return_index=True)

    return bp_dp, bp_x

def project_info(project=None, runid=None):
    outptu_dir = None

    if project == 'ICOGS':
        output_dir = f'/sciclone/schism10/feiye/ICOGS/{runid}/PostP/'
    elif project == 'STOFS3D-v4':
        output_dir = f'/sciclone/schism10/feiye/STOFS3D-v4/Outputs/O{runid.replace("RUN","")}/'
    else:
        print('Undefined project name')
        
    return output_dir

def get_datasets(project='STOFS3D-v4', runids=None):
    ths = []
    for runid in runids:
        runid1 = runid.split('RUN')[-1].strip()
        outdir = project_info(project=project, runid=runid)
        th = TimeHistory(
            file_name=f'{outdir}/elev.dat.missi.{runid}',
            start_time_str='2021-05-01 00:00:00', sec_per_time_unit=86400
        )
        ths.append(th)

    datasets = []
    for th, runid in zip(ths, runids):
        datasets.append(th.data[::1, :])
    
    return datasets


def transect_anim(gridfile=None, bpfile=None, project='STOFS3D-v4', runids=None):
    bp_dp, bp_x = get_transect_shape_from_grid(gridfile=gridfile, bpfile=bpfile)
    datasets = get_datasets(project=project, runids=runids)
    vis = testAnimation(datasets, runids, bp_x, bp_dp)
    return vis

def plot_snapshot(gridfile=None, bpfile=None, project=None, runids=[], snapshot_idx=-1):
    bp_dp, bp_x = get_transect_shape_from_grid(gridfile, bpfile)
    datasets = get_datasets(project='STOFS3D-v4', runids=runids)

    plt.plot(bp_x, -bp_dp/10)
    plt.plot(bp_x, datasets[0][0, :])
    for i, runid in enumerate(runids):
        plt.plot(bp_x, datasets[i][snapshot_idx, :], label=runid)
    plt.legend()
    plt.show()

# %%
if __name__ == "__main__":
    gridfile = '/sciclone/schism10/feiye/ICOGS/RUN10g/hgrid.npz'
    bpfile = '/sciclone/schism10/feiye/ICOGS/BPfiles/missi.bp'
    project = 'ICOGS'
    runids=['RUN10a', 'RUN10b', 'RUN10d', 'RUN10e', 'RUN10g']

    transect_anim(gridfile=gridfile, bpfile=bpfile, project=project, runids=runids).draw()

    rc('animation', html='jshtml')
    anim
