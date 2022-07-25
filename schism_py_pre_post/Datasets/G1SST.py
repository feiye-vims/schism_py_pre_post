from datetime import datetime, timedelta
import os
from matplotlib import pyplot
import numpy
from netCDF4 import Dataset as nc
import numpy as np


class G1SST():
    """ class for handling G1SST """
    def __init__(self, filename, basetime_str='1981-01-01T00:00:00Z'):
        self.source_file = filename
        self.nlon = 0
        self.nlat = 0
        self.lon = []
        self.lat = []
        self.depth = []
        self.time = []  # seconds since 1970-1-1
        self.t_datetime = []

        file_ext = filename.split(".")[-1]
        if file_ext == 'nc':
            if not os.path.exists(filename):
                if os.path.exists(filename + '.bz2'):
                    print(f'{filename} does not exist, extracting from {filename}.bz2')
                    os.system(f'bzip2 -dk {filename}.bz2')
                else:
                    raise Exception(f'{filename}* does not exist')
            file = nc(filename)
            t = file.variables['time']  # seconds since 1981
            base_time = datetime.strptime(t.units[-19:], "%Y-%m-%d %H:%M:%S")

            if len(t) > 1:
                t = numpy.squeeze(t)
                if t.ndim >= 2:
                    raise Exception("the time array has multiple dimensions")
            else:
                t = numpy.array(t)

            self.t_datetime = [base_time + timedelta(seconds=x.astype(float)) for x in t]

            start_time_str = self.t_datetime[0].strftime("%Y-%m-%d")
            end_time_str = self.t_datetime[-1].strftime("%Y-%m-%d")
            print("start_time: " + start_time_str)
            print("end_time: " + end_time_str)

            self.lon = numpy.array(file.variables['lon'])
            self.lat = numpy.array(file.variables['lat'])
        else:
            raise Exception('unkown file type')

    def get_var(self, var_str, shift=-273.15):
        file = nc(self.source_file)
        var = numpy.array(file.variables[var_str]) + shift
        return var

    def plot_init(
        self, this_var=None,
        t_idx=0, this_time_str=None,
        label_str=None, i_overwrite=True
    ):
        if this_time_str is not None:  # this_time_str prevails over t_idx
            t_idx = np.argmin(np.abs(np.array(self.t_datetime) - datetime.strptime(this_time_str, "%Y-%m-%d")))
        profile_time_str = datetime.strftime(self.t_datetime[t_idx], "%Y-%m-%d")

        if label_str is None:
            label_str = profile_time_str

        var = self.get_var(this_var)

        return [var, t_idx, profile_time_str, label_str]

    def plot_2d(
        self,
        i_show_plot=0,  # -1: don't draw; 0: draw but not show; 1: show
        this_var=None,
        t_idx=0, this_time_str=None,
        label_str=None, i_overwrite=True,
        xlim=None, ylim=None, clim=None,
        fig=None, model_basetime='2018-08-24T00:00:00Z'
    ):
        if xlim is None:
            xlim = [-83, -60]
        if ylim is None:
            ylim = [23, 43]
        if clim is None:
            clim = [16, 32]


        [var, t_idx, profile_time_str, label_str] =\
            self.plot_init(this_var, t_idx, this_time_str, label_str, i_overwrite)

        var_lyr = var[t_idx, :, :]
        var_lyr = numpy.where(var_lyr < -9999, numpy.nan, var_lyr)

        model_day_str = (self.t_datetime[t_idx] - datetime.strptime(model_basetime, "%Y-%m-%dT%H:%M:%SZ")).days

        pyplot.rcParams.update({'font.size': 15})
        if fig is None:
            fig, ax = pyplot.subplots(1, 1, figsize=(13, 8))

        dx = (self.lon[1]-self.lon[0])/2.
        dy = (self.lat[1]-self.lat[0])/2.
        extent = [self.lon[0]-dx, self.lon[-1]+dx, self.lat[0]-dy, self.lat[-1]+dy]
        c = ax.imshow(var_lyr, cmap='jet', vmin=clim[0], vmax=clim[1], extent=extent, interpolation="nearest", origin='lower')
        ax.title.set_text(f'{profile_time_str}; Day {model_day_str}')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        fig.colorbar(c, ax=ax)
        if i_show_plot == 0:
            pyplot.draw()
        if i_show_plot == 1:
            pyplot.show()

        return fig
