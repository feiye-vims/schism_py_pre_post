'''
Class for animating time series data
'''

from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from pylib import schism_grid
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, rc
from IPython.display import HTML

from schism_py_pre_post.Grid.Bpfile import Bpfile
import matplotlib

import plotly.io as pio


matplotlib.rcParams['animation.embed_limit'] = 2**128
pio.renderers.default = 'iframe_connected'

class testAnimation:
    def __init__(self, datasets, runids, bp_x, bp_dp, bp_levee_height):
        self.plot_x = bp_x
        self.bp_x = bp_x
        self.bp_dp = bp_dp / 5  # scaled
        self.bp_levee_height = bp_levee_height
        self.runids = runids
        self.bp_levee_height[self.bp_levee_height<-9999] = np.nan

        self.datasets = datasets
        self.nt = datasets[0].shape[0]

        # First set up the figure, the axis, and the plot element we want to animate
        self.fig, ax = plt.subplots()
        ax.set_facecolor('white')
        plt.close()
        ax.set_xlim(( min(bp_x), max(bp_x)))
        ax.set_ylim((-10, 12))
        self.fig.set_size_inches(12, 8)

        # set lines (elevation from different runs)
        markers = ['o', '*', '+' , 's', 'd']
        self.lines = []
        for i, [runid, marker] in enumerate(zip(runids, markers)):
            line, = ax.plot([], [], lw=1, marker=marker, markevery=5, label=runid)
            self.lines.append(line)

        # set other lines: bottom (river bed), levee height
        self.bottom, = ax.plot([], [], lw=2, label='river bed')
        self.levee, = ax.plot([], [], lw=2, label='levee')
        ax.legend()
        
        # set other decorations
        self.ax.legend()

    # initialization function: plot the background of each frame
    def init(self):
        for i, _ in enumerate(self.runids):
            self.lines[i].set_data([], [])
        self.bottom.set_data([], [])
        self.levee.set_data([], [])
        return (self.lines+[self.bottom]+[self.levee])

    # animation function. This is called sequentially  
    def animate(self, timestep):
        x = self.plot_x
        for k, data in enumerate(self.datasets):
            y = data[timestep, :]
            self.lines[k].set_data(x, y)
        self.bottom.set_data(self.bp_x, -self.bp_dp)
        self.levee.set_data(self.bp_x, -self.bp_levee_height)
        return (self.lines+[self.bottom]+[self.levee])

    def draw(self):
        global anim 
        anim = animation.FuncAnimation(self.fig, self.animate, init_func=self.init,
                                       frames=self.nt, interval=100, blit=True)