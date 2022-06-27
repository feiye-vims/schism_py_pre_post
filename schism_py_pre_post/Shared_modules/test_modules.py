#!/usr/bin/env python

"""
# module for pre and post processing
"""

from pylib import *
import csv
import re
import copy
from datetime import datetime, timedelta
import pickle
import statistics
import os
from pandas import read_csv, to_datetime, to_numeric
from matplotlib import pyplot
import numpy
import mplcursors  # 
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import math
from netCDF4 import Dataset as nc
from shutil import copy2
import subprocess
import gsw
import glob
import numpy as np

# from mpl_toolkits.basemap import Basemap


def pickle_save(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    """ Deserialize a file of pickled objects. """
    with open(filename, "rb") as f:
        data = pickle.load(f)
        return data
        # while True:
        #     try:
        #         yield pickle.load(f)
        #     except EOFError:
        #         break


def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))


def running_mean1(X, N):
    Y = X * numpy.nan
    if (X.ndim > 1):  # multi-dimensional
        for i, x in enumerate(X.transpose()):
            Y[:, i] = numpy.convolve(x, numpy.ones((N,))/N)[(N-1):]
    else:  # 1-d array
        Y = numpy.convolve(X, numpy.ones((N,))/N)[(N-1):]

    return Y


def running_mean(X, N):
    cumsum = numpy.cumsum(numpy.insert(X, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def timezone_dict(timezone_name=None, out_idx=0):
    if timezone_name is not None:
        tz_dict = {
            'CDT': ['America/Chicago', '-05:00'],
            'CST': ['America/Chicago', '-06:00'],
            'EDT': ['US/Eastern', '-04:00'],
            'EST': ['US/Eastern', '-05:00'],
            'MDT': ['US/Mountain', '-06:00'],
            'MST': ['US/Mountain', '-07:00'],
        }
        return tz_dict[timezone_name][out_idx]
    else:  # default to UTC
        if out_idx == 0:
            return 'UTC'
        else:
            return ' 00:00:00+00:00'


def schism_naming(var_name_str=None, schism_var_name=None):
    var_dict = {
        'temp': 'temperature',
        'hvel': 'horizontal velocity',
        'vert': 'vertical velocity',
        'zcor': 'z coordinate',
        'salt': 'salinity'
    }

    if var_name_str is None and schism_var_name is not None:
        return(var_dict[schism_var_name])

    found_key = []
    if var_name_str is not None and schism_var_name is None:
        for key in var_dict.keys():
            if var_name_str.lower() in key.lower():
                found_key.append(key)
        if len(found_key) == 0:
            raise Exception("Cannot find key")
        elif len(found_key) == 1:
            print(found_key)
        else:
            raise Exception("Multiple keys found" + found_key)

    return found_key[0]


#def datenum(d):
#    d = datetime.strptime(d, '%Y-%m-%dT%H:%M:%SZ')
#    return 366 + d.toordinal() + (d - datetime.fromordinal(d.toordinal())).total_seconds()/(24*60*60)
def datenum(d, date_fmt='%Y-%m-%dT%H:%M:%SZ'):
    d = datetime.strptime(d, date_fmt)
    return 366 + d.toordinal() + (d - datetime.fromordinal(d.toordinal())).total_seconds()/(24*60*60)

def datenum2datetime(t):
    return datetime.utcfromtimestamp((t - datenum('1970-1-1T00:00:00Z')) * 86400)

# def datestr(t):
#     return datetime.utcfromtimestamp((t - datenum('1970-1-1T00:00:00Z')) * 86400)
def datestr(t, date_fmt='%Y-%m-%dT%H:%M:%SZ'):
    return datenum2datetime(t).strftime(date_fmt)


def time_str_to_s_1970(time_str, time_str_fmt=None):
    epoch = datetime.strptime('1970-1-1T00:00:00Z', '%Y-%m-%dT%H:%M:%SZ').timestamp()
    # somehow 1 hour is missing, temporarily adding it manually
    if time_str_fmt is None:
        s_since_1970 = numpy.array(datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S').timestamp()) - epoch + 3600
    else:
        s_since_1970 = numpy.array(datetime.strptime(time_str, time_str_fmt).timestamp()) - epoch + 3600

    return s_since_1970


def extract_schism(bpfiles=None, var=None, rundir=None, i_comb=None, stacks=None, i_overwrite=True, i_all_level=False, ver=9):
    '''
    extract schism results based on one or more bpfiles
    ver: output version (e.g., ver 10 needs read_output10*)
    '''
    runid = os.path.basename(os.path.normpath(rundir))

    if i_all_level:
        script = f'/sciclone/home10/feiye/bin/read_output{ver}_xy_all_level.viz'
    else:
        script = f'/sciclone/home10/feiye/bin/read_output{ver}_xyz.viz'

    out_filenames = []

    if not os.path.exists(rundir):
        raise Exception("Schism dir not exist: " + rundir)
    if not os.path.exists(rundir + '/outputs/'):
        raise Exception("Schism outputs dir not exist: " + rundir)
    if not os.path.exists(rundir + '/extract_in/'):
        os.mkdir(rundir + "/extract_in/")
    if os.path.islink(rundir + '/PostP'):
        os.remove(rundir + '/PostP')
    if not os.path.exists(rundir + '/PostP/'):
        os.mkdir(rundir + "/PostP/")
    if not os.path.exists(rundir + '/outputs/vgrid.in'):
        copy2(rundir + '/vgrid.in', rundir + '/outputs/')

    for bpfile in bpfiles:
        bpname = os.path.basename(os.path.normpath(bpfile))
        bpname = bpname.split('.')[0]
        copy2(bpfile, rundir + '/outputs/station.bp')

        # make read_output*.in
        read_in = rundir + '/extract_in/' + var + '.in'
        with open(read_in, 'w') as fout:
            if ver==9:
                if i_comb:
                    fout.write("1\n2\n1\n")
                else:
                    fout.write("0\n1\n2.e9\n")
            fout.write("1\n1\n")
            fout.write(var + '\n')
            fout.write("1\n")
            fout.write(str(stacks[0]) + ' ' + str(stacks[1]) + '\n')
            fout.write("1\n")

        # call fortran script
        if i_all_level:
            out_filename = rundir + '/PostP/' + var + '.all_level.' + bpname + '.' + runid
        else:
            out_filename = rundir + '/PostP/' + var + '.dat.' + bpname + '.' + runid

        if var == 'hvel':
            out_filenames.append(out_filename.replace('hvel', 'u'))
            out_filenames.append(out_filename.replace('hvel', 'v'))
        else:
            out_filenames.append(out_filename)

        for this_file in out_filenames:
            if not os.path.exists(this_file) or i_overwrite:
                print(">>>> extracting schism results ...")
                p = subprocess.Popen(script, stdin=open(read_in), cwd=(rundir + '/outputs/'))
                p.wait()
                # input("Press Enter to continue...")
                if var == 'hvel':
                    os.system(f'mv {rundir}/outputs/fort.18 {out_filename.replace("hvel", "u")}')
                    os.system(f'mv {rundir}/outputs/fort.19 {out_filename.replace("hvel", "v")}')
                else:
                    os.system(f'mv {rundir}/outputs/fort.18 {out_filename}')
            else:
                print(">>>> schism results already extracted, skipping ...")

    return out_filenames

def extract_schism1(bpfiles=None, var=None, rundir=None, i_comb=None, stacks=None, i_overwrite=True, i_all_level=False, ver=9):
    '''
    extract schism results based on one or more bpfiles
    ver: output version (e.g., ver 10 needs read_output10*)
    '''
    runid = os.path.basename(os.path.normpath(rundir))

    if i_all_level:
        script = f'/sciclone/home10/feiye/bin/read_output{ver}_xy_all_level.viz'
    else:
        script = f'/sciclone/home10/feiye/bin/read_output{ver}_xyz.viz'

    out_filenames = []

    if not os.path.exists(rundir):
        raise Exception("Schism dir not exist: " + rundir)
    if not os.path.exists(rundir + '/outputs/'):
        raise Exception("Schism outputs dir not exist: " + rundir)
    if not os.path.exists(rundir + '/extract_in/'):
        os.mkdir(rundir + "/extract_in/")
    if os.path.islink(rundir + '/PostP'):
        os.remove(rundir + '/PostP')
    if not os.path.exists(rundir + '/PostP/'):
        os.mkdir(rundir + "/PostP/")
    if not os.path.exists(rundir + '/outputs/vgrid.in'):
        shutil.copy2(rundir + '/vgrid.in', rundir + '/outputs/')

    for bpfile in bpfiles:
        bpname = os.path.basename(os.path.normpath(bpfile))
        bpname = bpname.split('.')[0]
        shutil.copy2(bpfile, rundir + '/outputs/station.bp')

        # make read_output*.in
        read_in = rundir + '/extract_in/' + var + '.in'
        with open(read_in, 'w') as fout:
            if ver==9:
                if i_comb:
                    fout.write("1\n2\n10\n")
                else:
                    fout.write("0\n1\n2.e9\n")
            fout.write("1\n1\n")
            fout.write(var + '\n')
            fout.write("1\n")
            fout.write(str(stacks[0]) + ' ' + str(stacks[1]) + '\n')
            fout.write("1\n")

        # call fortran script
        if i_all_level:
            out_filename = rundir + '/PostP/' + var + '.all_level.' + bpname + '.' + runid
        else:
            out_filename = rundir + '/PostP/' + var + '.dat.' + bpname + '.' + runid

        if var == 'hvel':
            out_filenames.append(out_filename.replace('hvel', 'u'))
            out_filenames.append(out_filename.replace('hvel', 'v'))
        else:
            out_filenames.append(out_filename)

        for this_file in out_filenames:
            if not os.path.exists(this_file) or i_overwrite:
                print(">>>> extracting schism results ...")
                p = subprocess.Popen(script, stdin=open(read_in), cwd=(rundir + '/outputs/'))
                p.wait()
                # input("Press Enter to continue...")
                if var == 'hvel':
                    os.system(f'mv {rundir}/outputs/fort.18 {out_filename.replace("hvel", "u")}')
                    os.system(f'mv {rundir}/outputs/fort.19 {out_filename.replace("hvel", "v")}')
                else:
                    os.system(f'mv {rundir}/outputs/fort.18 {out_filename}')
            else:
                print(">>>> schism results already extracted, skipping ...")

    return out_filenames


class Argo():
    def __init__(self, filename=None, argo_repo_dir=None):
        if argo_repo_dir is None:
            argo_repo_dir = '/sciclone/data10/feiye/ARGO_REPO/'  # default download dir
        if isinstance(filename, str):
            if not os.path.exists(filename):
                raise Exception("File not exist: " + filename)
        elif isinstance(filename, list):
            [t_min_str, t_max_str, y_min, y_max, x_min, x_max] = filename
            filename = argo_repo_dir + ('_'.join(filename[:2])).replace(" ", "_") + "_" +\
                str(y_min) + "_" + str(y_max) + "_" +\
                str(x_min) + "_" + str(x_max) + '.nc'
            if not os.path.exists(filename):
                # sample url
                # http://www.ifremer.fr/erddap/tabledap/ArgoFloats.nc?
                # &time%3E=2017-08-04T00%3A00%3A00Z&time%3C=2017-09-14T14%3A48%3A20Z
                # &latitude%3E=25&latitude%3C=29&longitude%3E=-96&longitude%3C=-91
                url = "http://www.ifremer.fr/erddap/tabledap/ArgoFloats.nc?"\
                    + "&time%3E=" + t_min_str.replace(":", "%3A")\
                    + "&time%3C=" + t_max_str.replace(":", "%3A")\
                    + "&latitude%3E=" + str(y_min)\
                    + "&latitude%3C=" + str(y_max)\
                    + "&longitude%3E=" + str(x_min)\
                    + "&longitude%3C=" + str(x_max)
                print("Trying to download data for: " + "filename at " + url)
                os.system("curl '" + url + "' > " + filename)
        else:
            raise Exception("File name format wrong")

        self.source_file = filename


class Prop():
    """ class for *.prop or other similar formats"""
    def __init__(self, filename, number_of_groups):
        """Read a usgs data file and initialize from its content """
        self.source_file = filename
        self.n_group = number_of_groups
        self.np_group = []
        self.ip_group = []

        with open(self.source_file) as fin:
            for k in range(0, self.n_group):
                self.np_group.append(int(fin.readline().split()[0]))
                print("Points in Group " + str(k+1) + ": " + str(self.np_group))
                self.ip_group.append(numpy.empty((self.np_group[k]), dtype=int))
                for i in range(0, self.np_group[k]):
                    self.ip_group[k][i] = int(fin.readline())
                fin.readline()  # empty line between groups
                if self.np_group[k] > 0:
                    print("p first: " + str(self.ip_group[k][0]))
                    print("p last: " + str(self.ip_group[k][-1]))

    def get_max(self):
        """get max prop value"""

    def get_ele(self, pval):
        """return ele id with specified prop value"""


class Hgrid():
    # pylint: disable=too-many-instance-attributes
    """ class for manipulating hgrid.*"""

    def __init__(self, filename):
        """Read a hgrid file and initialize from its content """
        self.source_file = filename
        self.n_e = -1
        self.n_p = -1
        self.node_np = type(None)()
        self.x = type(None)()
        self.y = type(None)()
        self.dp = type(None)()
        self.lon = type(None)()
        self.lat = type(None)()
        self.ele_np = type(None)()
        self.ele = type(None)()
        self.i34 = type(None)()
        self.n_ob_total = type(None)()
        self.n_lb_total = type(None)()
        self.n_ob = type(None)()
        self.n_lb = type(None)()
        self.iob_np = type(None)()
        self.ilb_np = type(None)()
        self.got_nodes = False
        self.got_eles = False
        self.got_bnds = False
        self.ele_xyz = None
        self.alternate_file = f'{os.path.splitext(filename)[0]}.pickle'

        if os.path.exists(self.alternate_file):
            print(f'reading from alternative file: {self.alternate_file}')
            self.__dict__ = copy.deepcopy(pickle_load(self.alternate_file)).__dict__.copy()
            pass
        else:
            # read ne and np
            with open(filename, newline='') as fin:
                n_row = 0
                reader = csv.reader(fin, delimiter=' ', skipinitialspace=True)
                for row in reader:
                    n_row = n_row + 1
                    if n_row == 2:
                        [self.n_e, self.n_p] = list(map(int, row))
                        break
            print('ne: ' + str(self.n_e) + '; np: ' + str(self.n_p))

            if self.n_e <= 0 or self.n_p <= 0:
                raise Exception('ne/np <= 0: ' +
                                str(self.n_e) + ', ' + str(self.n_p))

    def save_pkl(self):
        pickle_save(self, self.alternate_file)

    def get_nodes(self):
        """ get nodes' x, y, z """
        self.node_np = numpy.empty((self.n_p, 4))
        self.x = numpy.empty((self.n_p))
        self.y = numpy.empty((self.n_p))
        self.y = numpy.empty((self.n_p))

        with open(self.source_file, newline='') as fin:
            n_row = 0
            n_p = 0
            reader = csv.reader(fin, delimiter=' ', skipinitialspace=True)
            for row in reader:
                n_row = n_row + 1
                if n_row > 2:
                    if n_row > 2+self.n_p:
                        break
                    self.node_np[n_p, :] = list(map(float, row))
                    n_p += 1

            self.x = self.node_np[:, 1]
            self.y = self.node_np[:, 2]
            self.dp = self.node_np[:, 3]

            print("First node:")
            print(int(self.node_np[0, 0]), self.node_np[0, 1:])
            print("Last node:")
            print(int(self.node_np[-1, 0]), self.node_np[-1, 1:])

        self.got_nodes = True

    def get_elements(self):
        """ get elements' nodes """
        self.ele_np = numpy.empty((self.n_e, 6), dtype=int)
        self.ele = numpy.empty((self.n_e, 4), dtype=int)
        self.i34 = numpy.empty((self.n_e), dtype=int)
        with open(self.source_file, newline='') as fin:
            n_row = 0
            n_e = 0
            reader = csv.reader(fin, delimiter=' ', skipinitialspace=True)
            for row in reader:
                n_row = n_row + 1
                if n_row > 2+self.n_p:
                    if n_row > 2+self.n_p+self.n_e:
                        break
                    tmp = list(map(int, row))
                    self.ele_np[n_e, 0:len(tmp)] = tmp
                    n_e += 1
            
            self.i34 = self.ele_np[:, 1]
            self.ele = self.ele_np[:, 2:6] - 1  # 0-based index

            print("First element:")
            print(self.ele_np[0, :2+self.ele_np[0, 1]])
            print("Last element:")
            print(self.ele_np[-1, :2+self.ele_np[0, 1]])

        self.got_eles = True

    def get_ele_centers(self):
        """ get coordinates of element centers"""
        if self.ele_xyz is not None:
            print('ele_centers already calculated')
            return

        element_list = range(0, self.n_e)
        element_list_np = numpy.array(element_list)

        if not self.got_nodes:
            self.get_nodes()

        if not self.got_eles:
            self.get_elements()

        print('element_list size: ' + str(element_list_np.size))

        self.ele_xyz = numpy.empty((self.n_e, 3), dtype=float)

        ie3 = self.ele_np[element_list, 1] == 3
        self.ele_xyz[ie3, 0] = numpy.mean(self.x[self.ele[ie3, :-1]], axis=1)
        self.ele_xyz[ie3, 1] = numpy.mean(self.y[self.ele[ie3, :-1]], axis=1)
        self.ele_xyz[ie3, 2] = numpy.mean(self.dp[self.ele[ie3, :-1]], axis=1)
        ie4 = self.ele_np[element_list, 1] == 4
        self.ele_xyz[ie4, 0] = numpy.mean(self.x[self.ele[ie4, :]], axis=1)
        self.ele_xyz[ie4, 1] = numpy.mean(self.y[self.ele[ie4, :]], axis=1)
        self.ele_xyz[ie4, 2] = numpy.mean(self.dp[self.ele[ie4, :]], axis=1)

        print("First element coor:")
        print(self.ele_xyz[0, :])
        print("Last element coor:")
        print(self.ele_xyz[-1, :])

        return self.ele_xyz

    def get_boundaries(self):
        """ get open and land bounaries
            land bounaries and island boundaries are not distinguished """
        self.n_ob = []
        self.n_lb = []

        with open(self.source_file) as fin:
            for i in range(0, 2+self.n_p+self.n_e):
                fin.readline()

            self.n_ob_total = int(fin.readline().split()[0])
            max_obn = int(fin.readline().split()[0])
            self.iob_np = numpy.empty((self.n_ob_total, max_obn), dtype=int)

            for i in range(0, self.n_ob_total):
                self.n_ob.append(int(fin.readline().split()[0]))
                for j in range(0, self.n_ob[i]):
                    self.iob_np[i, j] = int(fin.readline())

            self.n_lb_total = int(fin.readline().split()[0])
            max_nlb = int(fin.readline().split()[0])
            self.ilb_np = numpy.empty((self.n_lb_total, max_nlb), dtype=int)

            for i in range(0, self.n_lb_total):
                self.n_lb.append(int(fin.readline().split()[0]))
                for j in range(0, self.n_lb[i]):
                    self.ilb_np[i, j] = int(fin.readline())

            print("Open Boundary: %d" % (self.n_ob_total))
            print("Land/Island Boundary: %d" % (self.n_lb_total))

        self.got_bnds = True

    def find_node_from_xy(self, x, y, tol=None):
        """ find nearest from x,y """
        if not self.got_nodes:
            self.get_nodes()

        # default tol is for pin-pointing a node
        if tol is None:
            tol = numpy.sqrt(
                (self.node_np[1, 1]-self.node_np[2, 1])**2+(self.node_np[1, 2]-self.node_np[2, 2])**2
            ) * 0.01

        dist = numpy.sqrt((self.node_np[:, 1]-x)**2+(self.node_np[:, 2]-y)**2)
        if min(dist) > tol:
            raise Exception("Distance > tolerance")
        return numpy.argmin(dist)

    def find_ele_from_xy(self, x, y):
        """ find nearest from x,y """
        if self.ele_xyz is None:
            self.ele_xyz = self.get_ele_centers()
        dist = numpy.sqrt((self.ele_xyz[:, 0]-x)**2+(self.ele_xyz[:, 1]-y)**2)
        return numpy.argmin(dist)

    def draw_grid_boundary(self, fig=None, ishow=False):
        """ draw boundaries """
        if not self.got_nodes:
            self.get_nodes()

        if not self.got_bnds:
            self.get_boundaries()

        if fig is None:
            fig, ax = pyplot.subplots()
        else:
            ax = None

        out_xy = numpy.empty(shape=(0, 2), dtype=float)

        for i in range(0, self.n_ob_total):
            this_poly_x = self.node_np[self.iob_np[i, :self.n_ob[i]]-1, 1]
            this_poly_y = self.node_np[self.iob_np[i, :self.n_ob[i]]-1, 2]
            pyplot.plot(this_poly_x, this_poly_y, color='0.5')
            tmp = numpy.array([this_poly_x, this_poly_y])
            out_xy = numpy.vstack((out_xy, tmp.transpose()))

        for i in range(0, self.n_lb_total):
            this_poly_x = self.node_np[self.ilb_np[i, :self.n_lb[i]]-1, 1]
            this_poly_y = self.node_np[self.ilb_np[i, :self.n_lb[i]]-1, 2]
            pyplot.plot(this_poly_x, this_poly_y, color='0.5')
            if i == 0:
                tmp = numpy.array([this_poly_x, this_poly_y])
                out_xy = numpy.vstack((out_xy, tmp.transpose()))

        # pickle.dump(fig, open('FigureObject.pickle', 'wb'))
        if ishow:
            pyplot.show()

        return [fig, ax, out_xy]

    def in_grid_domain(self, coord_xy):
        """
        test if points inside grid initially used for the NOAA runs,
        which has one open bnd and one land bnd.
        Need to make it generic in the future
        """

        if not self.got_nodes:
            self.get_nodes()

        if not self.got_bnds:
            self.get_boundaries()

        poly_x = []
        poly_y = []
        for i in range(0, self.n_ob_total):
            poly_x.append(self.node_np[self.iob_np[i, :self.n_ob[i]]-1, 1])
            poly_y.append(self.node_np[self.iob_np[i, :self.n_ob[i]]-1, 2])

        for i in range(0, 1):
            poly_x.append(self.node_np[self.ilb_np[i, 1:self.n_lb[i]-1]-1, 1])
            poly_y.append(self.node_np[self.ilb_np[i, 1:self.n_lb[i]-1]-1, 2])

        poly_x = numpy.hstack(poly_x)
        poly_y = numpy.hstack(poly_y)
        poly = Polygon([[p_x, p_y] for p_x, p_y in zip(poly_x, poly_y)])

        i_pt_in_poly = []
        for this_point in coord_xy:
            point = Point(this_point[0], this_point[1])
            i_pt_in_poly.append(poly.contains(point))

        return i_pt_in_poly

    def writer(self, out_file_name):
        """ get nodes """
        if not self.got_nodes:
            self.get_nodes()
        if not self.got_eles:
            self.get_elements()
        # if not self.got_bnds:
        #     self.get_boundaries()

        """ write self """
        with open(out_file_name, 'w') as fout:
            fout.write("\n")
            fout.write(str(self.n_e) + ' ' + str(self.n_p) + '\n')
            for i, _ in enumerate(self.node_np):
                fout.write(str(int(self.node_np[i, 0])) + " " +
                           ' '.join(map(str, self.node_np[i, 1:])) +
                           "\n")
            for i, _ in enumerate(self.ele_np):
                fout.write(' '.join(map(str, self.ele_np[i, :2+self.ele_np[0, 1]])) + "\n")


class ObsDataSlot():

    __slots__ = ["stations", "fID2id", "saved_file"]

    def __init__(self, fname, from_raw_data=False):

        self.saved_file = fname
        self.fID2id = None

        if not (self.saved_file.endswith('pkl') or self.saved_file.endswith('npz')):
            raise Exception(f'file type not recognized: {os.path.basename(self.saved_file)}')

        self.stations = []

        dirname = os.path.dirname(self.saved_file)
        if not os.path.exists(dirname):
            print('reading existing database')
            raise Exception("Dir not exist: " + dirname)

        if os.path.exists(self.saved_file) and not from_raw_data:
            if self.saved_file.endswith('pkl'):
                self.stations, self.fID2id = pickle_load(self.saved_file).copy()
            elif self.saved_file.endswith('npz'):
                self.stations = loadz(self.saved_file).stations.tolist().copy()
                self.fID2id = loadz(self.saved_file).fID2id.tolist().copy()
        else:
            print('creating database from raw data')

            dirname = os.path.dirname(self.saved_file)
            raw_files = glob.glob(dirname + '/' + '*.compact')

            if len(raw_files) == 0:
                raise Exception("No raw files to work with under: " + dirname)

            for n, raw_file in enumerate(raw_files):
                print(f"----------------reading: {raw_file}------------------")
                this_station = Station(raw_file)
                this_station.get_lonlat()
                self.stations.append(this_station)
                # if n==10: break

            self.save(self.saved_file)

    def set_fID2id(self):
        self.fID2id = dict()
        for i, station in enumerate(self.stations):
            self.fID2id[station.id] = i

    def save(self):
        if self.saved_file.endswith('pkl'):
            pickle_save([self.stations, self.fID2id], self.saved_file)
        elif self.saved_file.endswith('npz'):
            S = npz_data()
            S.stations = self.stations
            S.fId2id = self.fId2id
            save_npz(self.saved_file, S)


class ObsData():

    def __init__(self, fname=None, from_raw_data=False):

        self.saved_file = fname
        self.fID2id = dict()
        self.stations = []

        if fname is None:
            return

        if not (self.saved_file.endswith('pkl') or self.saved_file.endswith('npz')):
            raise Exception(f'file type not recognized: {os.path.basename(self.saved_file)}')

        dirname = os.path.dirname(self.saved_file)
        if not os.path.exists(dirname):
            print('reading existing database')
            raise Exception("Dir not exist: " + dirname)

        if os.path.exists(self.saved_file) and not from_raw_data:
            if self.saved_file.endswith('pkl'):
                try:
                    self.stations, self.fID2id = pickle_load(self.saved_file).copy()
                except: 
                    self.stations = pickle_load(self.saved_file).copy()
                    self.set_fID2id()
            elif self.saved_file.endswith('npz'):
                self.stations = loadz(self.saved_file).stations.tolist().copy()
                self.fID2id = loadz(self.saved_file).fID2id.tolist().copy()
        else:
            print('creating database from raw data')

            dirname = os.path.dirname(self.saved_file)
            raw_files = glob.glob(dirname + '/' + '*.compact')

            if len(raw_files) == 0:
                raise Exception("No raw files to work with under: " + dirname)

            for n, raw_file in enumerate(raw_files):
                print(f"----------------reading: {raw_file}------------------")
                this_station = Station(raw_file)
                this_station.get_lonlat()
                self.stations.append(this_station)
                # if n==10: break

            self.set_fID2id()
            self.save(self.saved_file)

    def append(self, others=[]):
        if isinstance(others, list):
            for other in others:
                self.stations.append(other.stations)
        else:
            self.stations = self.stations + others.stations

        self.set_fID2id()

    def set_fID2id(self):
        self.fID2id = dict()
        for i, station in enumerate(self.stations):
            self.fID2id[station.id] = i

    def save(self, fname=None):
        self.saved_file = fname

        if self.saved_file.endswith('pkl'):
            pickle_save([self.stations, self.fID2id], self.saved_file)
        elif self.saved_file.endswith('npz'):
            S = npz_data()
            S.stations = self.stations
            S.fId2id = self.fId2id
            save_npz(self.saved_file, S)


class Station():

    def __init__(self, fname):
        import pandas as pd
        __slots__ = ["source_file", "id", "lon", "lat", "df"]

        self.source_file = fname
        self.id = os.path.basename(fname).split('_')[1]
        self.lon = -999
        self.lat = -999
        self.df = read_csv(fname)
        self.df.iloc[:, 1:] = self.df.iloc[:, 1:].apply(to_numeric, errors='coerce')
        self.df = self.df.set_index(pd.DatetimeIndex(pd.to_datetime(self.df["datetime"], utc=True)))
        # self.df["1994-01-01":"1994-02-01"].plot()
        # pyplot.show()
        pass

    def get_lonlat(self):
        import glob

        f_info = glob.glob(os.path.dirname(self.source_file) + '/' + '*info')
        if len(f_info) != 1:
            raise Exception('cannot determine info file to get lon/lat for station: ' + os.path.basename(self.source_file))
        else:
            df = read_csv(f_info[0], delim_whitespace=True,  index_col=False, dtype={'station_id':str})  # use the str of station id as index
            self.lon = df[df['station_id'] == self.id]['lon'].to_numpy()[0]
            self.lat = df[df['station_id'] == self.id]['lat'].to_numpy()[0]
            pass


class UsgsData():
    # pylint: disable-msg=too-many-locals
    """Read a usgs data file """
    def __init__(self, filename, timezone='EDT'):
        self.source_file = []
        self.station_id = []
        self.station_name = []
        self.var_dict = {}
        self.station_var = []
        self.header = []
        self.headerlines = 0
        self.timezone = timezone  # set default tz, will be overwritten if tz is specified in downloaded files

        # download file if not exist
        if isinstance(filename, str):
            if not os.path.exists(filename):
                raise Exception("File not exist: " + filename)
        elif isinstance(filename, list):
            [w_dir, data_type, station_id, var_id, begin_date, end_date] = filename
            filename = w_dir + '_'.join(filename[1:]) + '.rdb'
            if not os.path.exists(filename):
                url = "https://nwis.waterservices.usgs.gov/nwis/"\
                    + data_type + "/?"\
                    + "&format=rdb"\
                    + "&sites=" + station_id\
                    + "&startDT=" + begin_date\
                    + "&endDT=" + end_date\
                    + "&parameterCd=" + var_id
                print("Trying to download data for: '" + filename)
                os.system("curl '" + url + "' > " + filename)
        else:
            raise Exception("File name format wrong")

        self.source_file = filename

        # read metadata only
        with open(filename, newline='') as fin:
            n_row = 0
            n_sta = 0
            reader = csv.reader(fin, delimiter='\t')
            for row in reader:
                n_row = n_row + 1
                if re.search(r"^#\s+USGS\s+\d+", row[0]):
                    tmp = re.split(r'[;,\s]\s*', row[0])
                    self.station_id.append(tmp[2])
                    self.station_name.append('_'.join(tmp[3:]))
                elif re.search(r"^#\s+\d+\s+\d+", row[0]):
                    if re.search(r"^#\s+\d+\s+\d+\s+\d+", row[0]):
                        tmp = re.split(r'[;,\s]\s*', row[0])
                        print(tmp)
                        self.var_dict['_'.join(tmp[4:])] = '_'.join(tmp[1:4])
                    else:
                        tmp = re.split(r'[;,\s]\s*', row[0])
                        print(tmp)
                        self.var_dict['_'.join(tmp[3:])] = '_'.join(tmp[1:3])
                elif not re.search("^#", row[0]):  # field names
                    self.header = row
                    self.headerlines = n_row + 1
                    n_sta = n_sta + 1
                    break  # the rest is data

        if self.station_id == "":
            raise Exception("Station id not found in " + filename)

        print(self.station_id)
        print(self.station_name)
        print(self.var_dict)
        print(self.header)

    def locate_col_by_header(self, ident_list):
        """locate a data column by matching a list of strings"""
        n_found = 0
        found_entry = []
        found_key = []

        for idx, key in enumerate(self.header):
            i_found = 1
            for ident in ident_list:
                if "DO_NOT_USE".lower() in key.lower()\
                        or ident.lower() not in key.lower():  # case-insensitive match
                    i_found = 0
                    break
            if i_found == 1:
                n_found += 1
                found_key.append(key)
                found_entry.append(idx)
                print("Found key   [ " + str(key) +
                      ": " + self.header[idx] + " ]")

        if n_found == 0:
            raise Exception("None of the keys matches: " + ', '.join(ident_list[:]))
        if n_found > 1:
            print(found_key)
            print("Warning: Multiple matches")

        print("Index in Col: " + str(found_entry[0]))
        return [n_found, found_entry[0], found_key[0]]

    def locate_col(self, ident_list):
        """locate a data column by matching a list of strings"""
        n_found = 0
        found_entry = []
        found_key = []

        for key in self.var_dict:
            i_found = 1
            for ident in ident_list:
                ident_str = '_'.join(ident.split(' '))
                if "DO_NOT_USE".lower() in key.lower()\
                        or ident_str.lower() not in key.lower():  # case-insensitive match
                    i_found = 0
                    break

            if i_found == 1:
                n_found += 1
                found_key.append(key)
                found_entry.append(self.var_dict[key])
                print("Found key   [ " + str(key) +
                      ": " + self.var_dict[key] + " ]")

        if n_found == 0:
            print("Warning: None of the keys matches: " + ', '.join(ident_list[:]))
            return [n_found, -1, 'NAN']
        if n_found > 1:
            print(found_key)
            print("Warning: Multiple matches")

        i_col = []
        for this_entry in found_entry:
            i_col.append(self.header.index(this_entry))
        print("Index in Col: " + str(i_col))
        return [n_found, i_col, found_key]

    def get_field(self, ident_list):
        """ get a specified field """
        field = []
        # extract data col
        [_, i_data_col, _] = self.locate_col_by_header(ident_list)

        with open(self.source_file, 'r') as fin:
            reader = csv.reader(fin, delimiter='\t')
            n_row = 0
            for rows in reader:
                n_row += 1
                if n_row > self.headerlines:
                    field.append(rows[i_data_col])
                    # print(rows[i_data_col])

        return field

    def get_time_series(self, ident_list, iplot, unit_conv, i_demean=0, running_mean_window=0, shift=0.0):
        """ extract data """

        # extract time col
        if 'datetime' in self.header:
            i_time_col = self.header.index('datetime')
        else:
            return [-1, -1, -1, -1, -1, None]

        # extract time zone col
        i_tz_col = -1
        if 'tz_cd' in self.header:
            i_tz_col = self.header.index('tz_cd')
        else:
            print("<<<<<<Warning: time zone not found in " + self.source_file)
            # return [-1, -1, -1, -1, -1]

        # extract data col
        [n_found, i_data_col, data_name] = self.locate_col(ident_list)
        if n_found == 0:
            return [-1, -1, n_found, -1, -1, None]

        # save a copy
        out_file_name = self.source_file + "_" + '_'.join(ident_list) + ".compact"
        n_data_rec = 0
        with open(out_file_name, 'w') as fout:
            with open(self.source_file, 'r') as fin:
                reader = csv.reader(fin, delimiter='\t')
                n_row = 0
                fout.write("datetime, " + ", ".join(data_name).strip() + "\n")
                for rows in reader:
                    n_row += 1
                    if n_row > self.headerlines:
                        if n_data_rec == 0 and i_tz_col >= 0:  # get time zone
                            self.timezone = rows[i_tz_col].strip()
                        if i_tz_col == -1:
                            tz_str = ' 00:00:00+00:00'  # for dv
                        else:
                            tz_str = timezone_dict(rows[i_tz_col], 1)
                        n_data_rec += 1
                        fout.write(rows[i_time_col].strip() + tz_str + ', ' +
                                   ", ".join([rows[i] for i in i_data_col]) +
                                   "\n")
                        # time_str.append(rows[i_time_col])

        # deal with exceptions
        # os.system("sed -i 's/Dis/nan/g' " + out_file_name)
        # use "to_numeric, errors='coerce'" below

        # plot
        series = read_csv(out_file_name, header=0,
                          names=["time_stamp"]+data_name,
                          sep=r',')
        print(series.head())

        if iplot == 1:
            fig = pyplot.figure()

        dat_x = to_datetime(series['time_stamp'], utc=True)

        all_y = []
        for this_var in data_name:
            if series[this_var].dtype == 'O':
                series[this_var] = series[this_var].str.strip()
                series[this_var] = series[this_var].apply(to_numeric, errors='coerce')
            val_y = series[this_var].to_numpy().astype(numpy.float)
            all_y.append(val_y)

            if iplot == 1:
                label_str = "USGS_" + self.station_id[0] + "_" + '_'.join(ident_list)
                pyplot.plot(dat_x, val_y*unit_conv+shift, label=this_var.split('_')[0])
                fig.autofmt_xdate()

        # average
        # mean of multiple time series (e.g. bot/surf/mid) if they exist
        mean_y = numpy.nanmean(numpy.array(all_y), axis=0)
        if (i_demean == 1):
            # temporal mean
            mean_y1 = numpy.nanmean(mean_y, axis=0)
            mean_y -= mean_y1
        if (running_mean_window > 0):
            mean_y = running_mean(mean_y, running_mean_window)
            dat_x = dat_x[-mean_y.size:]
        first_y = numpy.array(all_y)[0][0]

        if iplot == 1:
            if n_found > 1:
                pyplot.plot(dat_x, mean_y*unit_conv+shift, label="average")
            pyplot.legend(loc="upper left")
            pyplot.show()
            pyplot.close(fig)

        if ('conductance' in data_name[0]):
            mean_y = gsw.conversions.SP_from_C(mean_y*unit_conv, 25, 0)
            first_y = gsw.conversions.SP_from_C(first_y*unit_conv, 25, 0)
        else:
            mean_y *= unit_conv
            first_y *= unit_conv

        if iplot == 1:
            numpy.savez((out_file_name + '.npz'),
                        method='plot',
                        args=(dat_x, mean_y+shift),
                        kwargs=dict(label=label_str))

        return [
            0, n_data_rec, n_found,
            numpy.nanmean(mean_y+shift), first_y,
            out_file_name + '.npz'
        ]

    def get_time_series_simple(self, ident, unit_conv):
        """ extract data """

        # extract time col
        if 'datetime' in self.header:
            i_time_col = self.header.index('datetime')
        else:
            return [-1, -1, -1, -1, -1, None]

        # extract time zone col
        i_tz_col = -1
        if 'tz_cd' in self.header:
            i_tz_col = self.header.index('tz_cd')
        else:
            print("<<<<<<Warning: time zone not found in " + self.source_file)

        # extract data col
        [n_found, i_data_col, data_name] = self.locate_col(ident)
        if n_found != 1:
            raise Exception(f'more than one column of data found for {ident}')

        # save a copy
        out_file_name = self.source_file + "_" + '_'.join(ident) + ".compact"
        n_data_rec = 0
        with open(out_file_name, 'w') as fout:
            with open(self.source_file, 'r') as fin:
                reader = csv.reader(fin, delimiter='\t')
                n_row = 0
                fout.write("datetime, " + ", ".join(data_name).strip() + "\n")
                for rows in reader:
                    n_row += 1
                    if n_row > self.headerlines:
                        if n_data_rec == 0 and i_tz_col >= 0:  # get time zone
                            self.timezone = rows[i_tz_col].strip()
                        if i_tz_col == -1:
                            tz_str = ' 00:00:00+00:00'  # for dv
                        else:
                            tz_str = timezone_dict(rows[i_tz_col], 1)
                        n_data_rec += 1
                        fout.write(rows[i_time_col].strip() + tz_str + ', ' +
                                   ", ".join([rows[i] for i in i_data_col]) +
                                   "\n")
                        # time_str.append(rows[i_time_col])

        # deal with exceptions
        # os.system("sed -i 's/Dis/nan/g' " + out_file_name)
        # use "to_numeric, errors='coerce'" below

        # plot
        df = read_csv(out_file_name, header=0,
                      names=["time_stamp"]+data_name,
                      sep=r',')
        df['time_stamp'] = to_datetime(df['time_stamp'], utc=True)
        df[data_name[0]] = df[data_name[0]] * unit_conv
        df.rename(columns={data_name[0]: f'{data_name[0]}*{unit_conv}'}, inplace=True)
        print(df.head())

        return df


class TimeHistory():

    from pandas import read_csv, to_datetime, to_numeric

    """Class for manipulating *.th"""

    def __init__(self, file_name, start_time_str, mask_val, sec_per_time_unit=1):
  
        if mask_val is None:
            mask_val = -99

        """initialize """
        self.source_file = file_name
        self.start_time_str = start_time_str
        self.n_time = -999
        self.n_station = -999
        self.time = None
        self.data = None
        self.delta_t = -999.9
        self.datetime = None

        self.df = read_csv(file_name, delim_whitespace=True, index_col=False, header=None)
        self.df.rename(columns={0: 'datetime'}, inplace=True)

        [self.n_time, self.n_station] = list(self.df.shape)
        self.n_station -= 1  # first col is time
        self.time = self.df['datetime'].to_numpy()
        self.datetime = [datetime.strptime(self.start_time_str, '%Y-%m-%d %H:%M:%S')
                         + timedelta(seconds=x*sec_per_time_unit) for x in self.time]
        self.df['datetime'] = self.datetime
        # mask invalid values

        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(mask_val)) < 1e-5] = numpy.nan
        print("Number of times: " + str(self.n_time))
        print("Number of stations: " + str(self.n_station))

        self.delta_t = self.time[1] - self.time[0]
        if self.delta_t < 50:  # probably in days
            self.time *= 86400.0
            self.delta_t *= 86400.0
            print("Time unit converted from days to seconds")
        else:
            print("Time unit: seconds")

    def get_running_mean(self, station_idx, window):
        """ sort time series by time average"""
        if station_idx is None:
            valid_idx = range(0, self.n_station)
        else:
            valid_idx = station_idx
        print(valid_idx)

        data_rm = running_mean(self.data[:, valid_idx], window)

        return [data_rm]

    def get_time_average(self, station_idx, start_time_str=None, end_time_str=None):
        if station_idx == []:
            valid_idx = range(0, self.n_station)
        else:
            valid_idx = station_idx
        print(valid_idx)

        if start_time_str is None:
            idx1 = 0
        else:
            idx1 = self.get_snapshot(start_time_str)

        if end_time_str is None:
            idx2 = self.n_time
        else:
            idx2 = self.get_snapshot(end_time_str)

        data_mean = numpy.mean(self.data[idx1:idx2, valid_idx], axis=0)

        return data_mean

    def sort_time_average(self, station_idx):
        """ sort time series by time average"""
        if station_idx == []:
            valid_idx = range(0, self.n_station)
        else:
            valid_idx = station_idx
        print(valid_idx)

        data_mean = numpy.mean(self.data[:, valid_idx], axis=0)
        id_sort = numpy.argsort(data_mean)
        print("Preview on sorted: ")
        print(data_mean[id_sort[-1:-11:-1]])

        return [id_sort, data_mean]

    def get_snapshot(self, this_time_str=None, this_time_s_1970=None):
        """ get a snapshot closest to the input time"""
        if this_time_s_1970 == "" or this_time_s_1970 is None:
            if this_time_str == "" or this_time_str is None:
                raise Exception("Failed to get time")
            this_time_s_1970 = time_str_to_s_1970(this_time_str)
        ts_start_s_1970 = time_str_to_s_1970(self.start_time_str)

        t_idx = numpy.argmin(numpy.abs(self.time +
                                       ts_start_s_1970 -
                                       this_time_s_1970))
        return(t_idx)

    def get_max_station(self):
        """ get max among all stations """
        data_sum = numpy.sum(self.data, axis=0)
        max_idx = numpy.argmax(data_sum)
        print("Station with the largest time-average value: " + str(max_idx))

        return max_idx
    
    def simple_plot(self, idx, **args):
        pyplot.plot(self.datetime, self.data[:, idx], **args)
        

    def plot_ts(self, idx, title_str, i_plot=0, i_demean=0, running_mean_window=0, shift=0.0):
        """ plot specific col indices """
        """ idx is for data, not counting the first column of times"""
        if idx is None:
            idx = 0

        iMultiLine = False
        if isinstance(idx, (list, numpy.ndarray)):
            if (len(idx) > 1):
                iMultiLine = True
            else:
                idx = idx[0]

        fig, ax = pyplot.subplots()
        label_str = '_'.join(title_str.split('_')[0:2])

        # save a copy before further processing
        data_copy = self.data[:, idx]
        if i_demean == 1:
            self.data[:, idx] -= numpy.nanmean(data_copy, 0)
        elif i_demean == -1:  # hack, use a constant (mean)
            self.data[:, idx] = numpy.nanmean(data_copy, 0)

        if running_mean_window > 0:
            data_rm = running_mean(data_copy, running_mean_window)
            self.data[-data_rm.size:, idx] = data_rm

        if self.start_time_str == "":
            ax.plot(self.time, self.data[:, idx]+shift, label=label_str)
        else:
            time_origin_utc = datetime.strptime(self.start_time_str, '%Y-%m-%d %H:%M:%S')
            time_with_origin = self.time*10**6 + \
                self.timestamp_microsecond(datetime(1970, 1, 1), time_origin_utc)
            time_with_origin = time_with_origin/10**6  # to seconds
            time_with_origin_str = to_datetime(time_with_origin, unit='s')

            print("Start: ")
            print(time_with_origin_str[0])
            print("End: ")
            print(time_with_origin_str[-1])
            fig.autofmt_xdate()
            ax.plot(time_with_origin_str, self.data[:, idx]+shift, label=label_str)

        ax.legend(loc="upper left")

        if not iMultiLine:
            pickle.dump(
                fig, open(self.source_file + '_' + str(idx) + '_' + title_str + '.pickle', 'wb')
            )

        mplcursors.cursor(fig)
        if i_plot == 1:
            pyplot.show()
        pyplot.close(fig)

        if not iMultiLine:
            out_file_name = (self.source_file + '_' + title_str + '.npz')
            numpy.savez(out_file_name,
                        method='plot',
                        args=(time_with_origin_str, self.data[:, idx]+shift),
                        kwargs=dict(label=label_str))

        # restore data
        self.data[:, idx] = data_copy
        return out_file_name

    def time_unit_conversion(self, unit_conv):
        """convert time unit"""
        """Col 0 is time: unit_conv from sec to day is 1/86400"""
        self.time = self.time * unit_conv

    def data_unit_conversion(self, unit_conv):
        """convert data unit"""
        self.data = self.data * unit_conv

    def export_subset(self, station_idx, time_idx, i_reset_time, subset_filename):
        """extract a subset from the original *.th"""
        """by station_idx and time_idx"""
        """[] idx means no operation"""
        """if i_reset_time == 1, then the subset *.th starts from 0 time"""
        station_idx = numpy.array(station_idx)
        time_idx = numpy.array(time_idx)

        if (len(station_idx) > 0):
            if (len(time_idx) > 0):
                sub_data = self.data[time_idx, station_idx]
            else:
                sub_data = self.data[:, station_idx]

        if (len(time_idx) > 0):
            sub_time = self.time[time_idx]
        else:
            sub_time = self.time

        if (i_reset_time > 0):
            sub_time = sub_time - sub_time[0]

        with open(subset_filename, 'w') as fout:
            for i, _ in enumerate(sub_time):
                fout.write(str(sub_time[i]) + " " +
                           ' '.join(map(str, sub_data[i, :])) +
                           "\n")

    def writer(self, out_file_name):
        """ write self """
        with open(out_file_name, 'w') as fout:
            for i, _ in enumerate(self.time):
                fout.write(str(self.time[i]) + " " +
                           ' '.join(map(str, self.data[i, :])) +
                           "\n")

    @classmethod
    def timestamp_microsecond(cls, epoch, utc_time):
        """ convert date to microseconds """
        time_diff = utc_time - epoch
        assert time_diff.resolution == timedelta(microseconds=1)
        return (time_diff.days * 86400 + time_diff.seconds) * 10**6 + time_diff.microseconds


class Misc():  # pylint: disable=too-few-public-methods
    """class for misc functions"""


class GliderData():
    """class for glider data"""
    def __init__(self, filename, xyzt_str):
        self.source_file = filename
        self.n_data_points = 0
        self.start_time_str = []
        self.end_time_str = []
        self.t_datetime = []
        self.xyzt = []
        self.n_profiles = 0
        self.profile_idx = []
        self.profile_ids = []
        self.profile_ids_by_point = []
        self.profile_xy = []
        self.profile_time = []
        self.sorted_profile_time_idx = []
        self.i2d = False

        file_ext = filename.split(".")[-1]
        if file_ext == 'nc':
            file = nc(filename)
            t = file.variables[xyzt_str[3]][:]  # seconds since 1970-1-1
            t = numpy.squeeze(t)

            if t.ndim >= 2:
                self.i2d = True
                [self.n_profiles, n_depth] = t.shape
                t_2d_mask = numpy.ma.getmaskarray(t)
                t = t.flatten()
                t_mask = numpy.ma.getmaskarray(t)
                t = t[~t_mask]

            self.t_datetime = to_datetime(t, unit='s')
            start_time_str = to_datetime(t[0], unit='s')
            end_time_str = to_datetime(t[-1], unit='s')
            print("start_time: " + str(start_time_str))
            print("end_time: " + str(end_time_str))

            self.xyzt = numpy.empty((len(t), 4))
            self.xyzt[:, 3] = t  # seconds since 1970-1-1
            if not self.i2d:
                self.xyzt[:, 0] = file.variables[xyzt_str[0]][:]
                self.xyzt[:, 1] = file.variables[xyzt_str[1]][:]
                depth = file.variables[xyzt_str[2]][:]
                d_mask = numpy.ma.getmaskarray(depth)
                depth[d_mask] = numpy.nan
                self.xyzt[:, 2] = depth

                if xyzt_str[4] != '':
                    n_points_in_profiles = file.variables[xyzt_str[4]][:]
                    self.n_profiles = len(n_points_in_profiles)
                    self.profile_idx = numpy.empty((self.n_profiles), dtype=int)
                    self.profile_idx[0] = 0
                    for i in range(1, self.n_profiles):
                        self.profile_idx[i] = self.profile_idx[i-1]+n_points_in_profiles[i-1]

                if xyzt_str[5] != '':
                    if len(xyzt_str[5]) == 1:
                        self.profile_ids = file.variables[xyzt_str[5][0]][:]
                    elif len(xyzt_str[5]) == 2:
                        tmp1 = file.variables[xyzt_str[5][0]][:]
                        tmp2 = file.variables[xyzt_str[5][1]][:]
                        for str1, str2 in zip(tmp1, tmp2):
                            self.profile_ids.append(str1 + "_" + str(str2))
                    else:
                        raise Exception('failed to parse profile_ids')
                    if len(self.profile_idx) == 0:
                        if len(self.profile_ids) == len(depth):
                            self.profile_ids_by_point = self.profile_ids
                            self.profile_ids = numpy.unique(self.profile_ids_by_point)
                            self.n_profiles = len(self.profile_ids)
                            self.profile_idx = numpy.empty((self.n_profiles), dtype=int)

                            pid_last = -999
                            n = 0
                            for i, pid in enumerate(self.profile_ids_by_point):
                                if pid != pid_last:
                                    self.profile_idx[n] = i
                                    n += 1
                                    pid_last = pid
                            if n != self.n_profiles:
                                raise Exception("Profile miscount")
                        else:
                            raise Exception("failed to get profile_idx from profile_ids")

                if len(self.profile_idx) == 0 or len(self.profile_ids) == 0:
                    raise Exception("failed to get profile_ids")

            else:
                self.xyzt[:, 0] = file.variables[xyzt_str[0]][:].flatten()[~t_mask]
                self.xyzt[:, 1] = file.variables[xyzt_str[1]][:].flatten()[~t_mask]
                depth = file.variables[xyzt_str[2]][:].flatten()[~t_mask]
                # depth[depth <= -9999] = numpy.nan
                self.xyzt[:, 2] = depth

                # get profile_idx
                self.profile_idx = numpy.empty((self.n_profiles), dtype=int)
                self.profile_idx[0] = 0
                for i in range(1, self.n_profiles):
                    this_mask = t_2d_mask[i, :]
                    self.profile_idx[i] = self.profile_idx[i-1] + len(this_mask[~this_mask])
                    if len(this_mask[~this_mask]) == 0:
                        print("Warning: empty profile: " + str(i) + "; " + str(self.profile_ids[i]))

                if xyzt_str[5] != '':
                    self.profile_ids = file.variables[xyzt_str[5]][:].flatten()
                else:
                    raise Exception("failed to get profile_ids")

            # assign profile id on each point
            if len(self.profile_ids_by_point) == 0:
                self.profile_ids_by_point = numpy.empty((len(depth)), dtype=int)
                for i in range(0, self.n_profiles-1):
                    idx1 = self.profile_idx[i]
                    idx2 = self.profile_idx[i+1]
                    self.profile_ids_by_point[idx1:idx2] = self.profile_ids[i]
                self.profile_ids_by_point[self.profile_idx[-1]:] = self.profile_ids[-1]

        elif file_ext == 'csv':
            with open(filename, newline='') as fin:
                df = read_csv(fin, sep=',')

            if xyzt_str[5] != '':
                self.profile_ids_by_point = numpy.array(df[xyzt_str[5]].to_numpy()[1:], dtype=int)

                pid_last = -999
                self.n_profiles = 0
                for i, pid in enumerate(self.profile_ids_by_point):
                    if pid != pid_last:
                        self.n_profiles += 1
                        pid_last = pid
                self.profile_idx = numpy.empty(self.n_profiles, dtype=int)

                pid_last = -999
                n = 0
                for i, pid in enumerate(self.profile_ids_by_point):
                    if pid != pid_last:
                        self.profile_idx[n] = i
                        n += 1
                        pid_last = pid
                if n != self.n_profiles:
                    raise Exception("Profile miscount")
            else:
                raise Exception("Profile string undefined")
            t = df[xyzt_str[3]].to_numpy()[1:]
            epoch = datetime.strptime('1970-1-1T00:00:00Z', '%Y-%m-%dT%H:%M:%SZ').timestamp()
            # somehow 1 hour is missing, temporarily adding it manually
            t_s_since_1970 = numpy.array([datetime.strptime(this_t, '%Y-%m-%dT%H:%M:%SZ')
                                         .timestamp() for this_t in t]) - epoch + 3600
            self.t_datetime = to_datetime(t_s_since_1970, unit='s')
            self.xyzt = numpy.empty((len(t), 4))
            self.xyzt[:, 3] = t_s_since_1970
            for i in range(0, 3):
                self.xyzt[:, i] = df[xyzt_str[i]].to_numpy()[1:]
        else:
            raise Exception("unknown file type for glider data")

        # get some additional info for each profile
        self.profile_time = numpy.empty((self.n_profiles), dtype=float)
        self.profile_xy = numpy.empty((self.n_profiles, 2), dtype=float)
        for i in range(0, self.n_profiles-1):
            idx1 = self.profile_idx[i]
            idx2 = self.profile_idx[i+1]
            self.profile_xy[i, 0] = numpy.mean(self.xyzt[idx1:idx2, 0])
            self.profile_xy[i, 1] = numpy.mean(self.xyzt[idx1:idx2, 1])
            self.profile_time[i] = numpy.mean(self.xyzt[idx1:idx2, 3])
        self.profile_xy[-1, 0] = numpy.mean(self.xyzt[self.profile_idx[-1]:, 0])
        self.profile_xy[-1, 1] = numpy.mean(self.xyzt[self.profile_idx[-1]:, 1])
        self.profile_time[-1] = numpy.mean(self.xyzt[self.profile_idx[-1]:, 3])

        self.sorted_profile_time_idx = self.profile_time.argsort()

    def get_var(self, var_name):
        """ get a specific variable """
        file = nc(self.source_file)
        this_var = file.variables[var_name][:]
        print("Got " + var_name + " from " + self.source_file)

        mask = numpy.ma.getmaskarray(this_var)
        this_var[mask] = numpy.nan

        return numpy.array(this_var)

    def get_single_profile(self, idx, varname):
        start_idx = self.profile_idx[idx]
        if idx == self.n_profiles-1:
            end_idx = self.n_data_points
        else:
            end_idx = self.profile_idx[idx+1]

        lon = self.xyzt[start_idx:end_idx, 0]
        lat = self.xyzt[start_idx:end_idx, 1]
        depth = self.xyzt[start_idx:end_idx, 2]
        time_str = self.t_datetime[start_idx:end_idx]
        time = self.xyzt[start_idx:end_idx, 3]

        tmp = self.get_var(varname)
        var = tmp[start_idx:end_idx]

        return[lon, lat, depth, time, time_str, var]

    def plot_single_profile(self, idx, i_show_plot, varname):
        """ plot single profile, idx starting from 0 """
        profile_time_str = str(to_datetime(self.profile_time[idx], unit='s'))
        out_file_name = (self.source_file + "."
                         + varname + "."
                         + profile_time_str.replace(" ", "_")
                         + ".npz")

        fig, ax = pyplot.subplots(1, 1, figsize=(13, 8))

        [_, _, depth, _, time_str, var] = self.get_single_profile(idx, varname)
        if numpy.isnan(depth).all():
            i_valid = False
        else:
            i_valid = True
        if varname == "":  # plot time against time if no var specified
            dat_x = time_str
            dat_y = -depth
        else:  # plot var against depth
            dat_x = var
            dat_y = -depth

        numpy.savez((out_file_name),
                    method='plot',
                    args=(dat_x, dat_y),
                    kwargs=dict(marker='o', label=profile_time_str))
                    # kwargs=dict(label=profile_time_str))

        # ax.plot(dat_x, dat_y, marker='o')
        ax.plot(dat_x, dat_y)
        ax.title.set_text(profile_time_str)
        if i_show_plot == 1:
            pyplot.show(fig)
        pyplot.close(fig)

        return [fig, ax, out_file_name, profile_time_str, i_valid]

    def get_profile_from_time(self, timestamp):
        """ get the id (base 0) on the closest timestamp"""
        epoch = datetime.strptime('1970-1-1T00:00:00Z', '%Y-%m-%dT%H:%M:%SZ').timestamp()
        # somehow 1 hour is missing, temporarily adding it manually
        t_s_since_1970 = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S').timestamp()-epoch+3600

        # point_idx = numpy.argmin(numpy.abs(self.xyzt[:, 3]-t_s_since_1970))
        # profile_id = self.profile_ids_by_point[point_idx]
        # profile_idx = numpy.where(self.profile_ids == profile_id)
        #
        # if len(profile_idx) != 1:
        #     raise Exception("failed to find a single profile index")
        #
        # print("Profile " + str(profile_id) + " at Index " + str(profile_idx[0]))
        # return [int(profile_idx[0]), profile_id]

        profile_idx = numpy.argmin(numpy.abs(self.profile_time-t_s_since_1970))
        profile_id = self.profile_ids[profile_idx]
        return [int(profile_idx), profile_id]

    @classmethod
    def timestamp_microsecond(cls, epoch, utc_time):
        """ convert date to microseconds """
        time_diff = utc_time - epoch
        assert time_diff.resolution == timedelta(microseconds=1)
        return (time_diff.days * 86400 + time_diff.seconds) * 10**6 + time_diff.microseconds


class Bpfile():

    """class for bpfiles"""

    def __init__(self, filename, cols=4):
        """Initialization """
        self.nodes = []
        self.n_nodes = 0
        self.info = []
        self.ncol = cols

        self.reader(filename, cols)
        self.make_dataframe()

        print(str(self.n_nodes) + "\n")
        print(str(self.nodes[0][0]) + " " + str(self.nodes[0][1]) + "\n")

    def reader(self, file_name, cols):
        """ Read existing *.reg """
        with open(file_name) as fin:
            fin.readline()
            self.n_nodes = int(fin.readline().split()[0])
            for _ in range(0, self.n_nodes):
                line = fin.readline()
                self.nodes.append(list(map(float, line.split()[slice(1, 4)])))
                self.info.append(list(line.split()[slice(4, cols)]))
        self.nodes = numpy.array(self.nodes, dtype=float)

    def writer(self, out_file_name, ncol=2):
        """ nodes should be shaped as [:n_nodes][:2] """
        with open(out_file_name, 'w') as fout:
            fout.write("\n")
            fout.write(str(len(self.nodes)) + "\n")
            for i, _ in enumerate(self.nodes):
                fout.write(str(i+1) + " ")
                for j in range(0, ncol):
                    fout.write(str(self.nodes[i][j]) + " ")
                fout.write("\n")
    
    def make_dataframe(self, row_name=['lon', 'lat', 'z'], col_name=None):
        import pandas as pd
        if col_name is None:
            self.st_id = []
            for i, st in enumerate(self.info):
                if self.ncol < 4:
                    self.st_id.append(f'{i}')
                else:
                    if len(st)>0:  # station name
                        self.st_id.append(st[0][1:])
                    else:
                        self.st_id.append(f'Station_{i+1}')
            col_name = self.st_id

        self.df = pd.DataFrame(data=self.nodes.T, index=row_name, columns=col_name)
        return self.df


class XmRegion():

    """class for regions in xmgredit5"""

    def __init__(self, input_region):
        """Initialization """
        self.nodes = []
        self.n_nodes = 0

        if isinstance(input_region, str):  # Read existing *.reg
            self.reader(input_region)
        elif isinstance(input_region, (list, numpy.ndarray)):
            if len(input_region[0]) != 2:
                raise Exception("invalid input region list")
            self.nodes = input_region
            self.n_nodes = len(input_region)
        else:
            print(type(input_region))
            raise Exception("invalid input region type")

        print(str(self.n_nodes) + "\n")
        print(str(self.nodes[0][0]) + " " + str(self.nodes[0][1]) + "\n")

    def reader(self, file_name):
        """ Read existing *.reg """
        with open(file_name) as fin:
            fin.readline()
            fin.readline()
            self.n_nodes = int(fin.readline().split()[0])
            for _ in range(0, self.n_nodes):
                self.nodes.append(list(map(float, fin.readline().split()[0:2])))

    def writer(self, out_file_name):
        """ nodes should be shaped as [:n_nodes][:2] """
        with open(out_file_name, 'w') as fout:
            fout.write("\n")
            fout.write("1\n")
            fout.write(str(len(self.nodes)) + " 1 \n")
            for i, _ in enumerate(self.nodes):
                fout.write(str(self.nodes[i][0]) + " " + str(self.nodes[i][1]) + "\n")


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


class HYCOM():
    """ class for handling HYCOM"""

    def __init__(self, filename):
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
            file = nc(filename)
            t = file.variables['time']  # hours since 2000-1-1
            t = numpy.squeeze(t)

            if t.ndim >= 2:
                raise Exception("the time array has multiple dimensions")

            base_time_s_1970 = time_str_to_s_1970('2000-01-01 00:00:00')
            self.time = t*3600 + base_time_s_1970  # convert to seconds then add shift
            self.t_datetime = to_datetime(self.time, unit='s')
            start_time_str = to_datetime(self.time[0], unit='s')
            end_time_str = to_datetime(self.time[-1], unit='s')
            print("start_time: " + str(start_time_str))
            print("end_time: " + str(end_time_str))

            try:
                self.lon = numpy.array(file.variables['xlon'])
            except KeyError:
                self.lon = numpy.array(file.variables['lon'])-360.0
            try:
                self.lat = numpy.array(file.variables['ylat'])
            except KeyError:
                self.lat = numpy.array(file.variables['lat'])
            try:
                self.depth = numpy.array(file.variables['depth'])
            except Exception:
                self.depth = None
        else:
            raise Exception('unkown file type')

    def get_time(self, this_time_str=None, this_time_s_1970=None):
        if this_time_s_1970 is None:
            if this_time_str is None:
                raise Exception("Failed to get time")
            this_time_s_1970 = time_str_to_s_1970(this_time_str)

        idx = numpy.argmin(numpy.abs(self.time - this_time_s_1970))

        rat = [0.0, 0.0]
        indices = [-1, -1]
        if self.time[idx] <= this_time_s_1970:
            rat[1] = (this_time_s_1970 - self.time[idx])/(self.time[idx+1] - self.time[idx])
            rat[0] = 1.0 - rat[1]
            indices = [idx, idx+1]
        else:
            rat[0] = (this_time_s_1970 - self.time[idx-1])/(self.time[idx] - self.time[idx-1])
            rat[1] = 1.0 - rat[0]
            indices = [idx-1, idx]

        return [rat, indices]

    def get_depth(self, this_depth):
        idx = numpy.argmin(numpy.abs(self.depth - this_depth))

        rat = [0.0, 0.0]
        indices = [-1, -1]
        if self.depth[idx] <= this_depth:
            rat[1] = (this_depth - self.depth[idx])/(self.depth[idx+1] - self.depth[idx])
            rat[0] = 1.0 - rat[1]
            indices = [idx, idx+1]
        else:
            rat[0] = (this_depth - self.depth[idx-1])/(self.depth[idx] - self.depth[idx-1])
            rat[1] = 1.0 - rat[0]
            indices = [idx-1, idx]

        return [rat, indices]

    def get_var(self, var_str):
        file = nc(self.source_file)
        var = numpy.array(file.variables[var_str])
        return var

    def plot_init(
        self, this_var=None,
        this_time_str=None, this_time_s_1970=None,
        label_str=None, i_overwrite=True
    ):
        if this_time_s_1970 is None:
            if this_time_str is None:
                raise Exception("Failed to get time")
            this_time_s_1970 = time_str_to_s_1970(this_time_str)

        t_idx = numpy.argmin(numpy.abs(self.time - this_time_s_1970))
        profile_time_str = str(self.t_datetime[t_idx])

        if label_str is None:
            label_str = profile_time_str

        out_file_name = (self.source_file + "."
                         + this_var + "."
                         + profile_time_str.replace(" ", "_")
                         + ".npz")
        if not os.path.exists(out_file_name) or i_overwrite:
            # reading file contents
            var = self.get_var(this_var)
            i_exist = False
        else:
            var = None
            i_exist = True

        return [var, t_idx, profile_time_str, out_file_name, label_str, i_exist]

    def plot_2d(
        self,
        i_show_plot=0,  # -1: don't draw; 0: draw but not show; 1: show
        this_var=None,
        this_time_str=None, this_time_s_1970=None,
        i_overwrite=True,
        xlim=None, ylim=None, clim=None,
        lyr=0,
        label_str=None,
        fig=None
    ):
        if xlim is None:
            xlim = [-90, -62]
        if ylim is None:
            ylim = [20, 42]
        if clim is None:
            clim = [26, 40]

        [var, t_idx, profile_time_str, out_file_name, label_str, i_exist] =\
            self.plot_init(this_var, this_time_str,
                           this_time_s_1970, label_str, i_overwrite)

        # skip existing plot and return now
        if (not i_overwrite) and (i_exist) and (i_show_plot == -1):
            return out_file_name

        var_lyr = var[t_idx, lyr, :, :]
        var_lyr = numpy.where(var_lyr < -9999, numpy.nan, var_lyr)

        if fig is None:
            fig, ax = pyplot.subplots(1, 1, figsize=(13, 8))

        c = ax.pcolor(self.lon, self.lat, var_lyr, cmap='jet', vmin=clim[0], vmax=clim[1])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        fig.colorbar(c, ax=ax)
        if i_show_plot == 0:
            pyplot.draw()
        if i_show_plot == 1:
            pyplot.show(fig)

        return fig

    def plot_profile(
        self, i_show_plot=False, this_var=None,
        this_lon=None, this_lat=None,
        this_time_str=None, this_time_s_1970=None,
        label_str=None,
        i_overwrite=False
    ):
        """ get a vertical profile of a specified variable at the specified location and time"""
        [var, t_idx, profile_time_str, out_file_name, label_str, i_exist] =\
            self.plot_init(this_var, this_time_str,
                           this_time_s_1970, label_str, i_overwrite)

        # skip existing plot and return now
        if (not i_overwrite) and (i_exist) and (not i_show_plot):
            return out_file_name

        lon_idx = numpy.argmin(numpy.abs(self.lon - this_lon))
        lat_idx = numpy.argmin(numpy.abs(self.lat - this_lat))

        fig, ax = pyplot.subplots(1, 1, figsize=(13, 8))
        dat_x = var[t_idx, :, lat_idx, lon_idx]
        dat_x = numpy.where(dat_x < -9999, numpy.nan, dat_x)
        dat_y = -self.depth

        numpy.savez((out_file_name),
                    method='plot',
                    args=(dat_x, dat_y),
                    kwargs=dict(label=label_str))
        #           kwargs=dict(marker='o', label=profile_time_str))

        # ax.plot(dat_x, dat_y, marker='o')
        ax.plot(dat_x, dat_y)
        ax.title.set_text(profile_time_str)
        if i_show_plot:
            pyplot.show(fig)
        pyplot.close(fig)

        return out_file_name


class CBP_casts():

    """ class for handling CBP casts"""

    def __init__(self, filename=None, data=None):
        if filename is not None:
            self.source_file = filename
            self.data = read_csv(self.source_file)
            my_datetime = to_datetime(
                self.data['SampleDate'] + " " + self.data['SampleTime'], format='%m/%d/%Y %H:%M:%S'
            ).dt.tz_localize(timezone_dict('EDT')).dt.tz_convert('GMT').dt.tz_convert(None)
            self.data['Datetime'] = my_datetime
        elif data is not None:
            self.source_file = None
            self.data = data

    def get_station(self, station_name):
        return CBP_casts(None, self.data.query('Station==@station_name'))

    def get_var(self, var_name):
        return CBP_casts(None, self.data.query('Parameter==@var_name'))

    def get_snapshot(self, datetime_str, datetime_str_fmt=None):
        if datetime_str_fmt is None:
            nearest_time = nearest(self.data.Datetime, to_datetime(datetime_str))
        else:
            nearest_time = nearest(self.datetime, to_datetime(datetime_str, format=datetime_str_fmt))
        print('nearest time: ' + nearest_time.strftime('%Y-%m-%d %H:%M:%S'))
        return CBP_casts(None, self.data.query('Datetime==@nearest_time'))

    def sort_by(self, var_str):
        if var_str in self.data.columns:
            self.data = self.data.sort_values(by=[var_str])
            self.data = self.data.drop_duplicates(subset=var_str, keep="first")
            print('CBP data sorted by ' + var_str)
            return CBP_casts(None, self.data)
        else:
            raise Exception(var_str + ' is not in CBP keys')

    def write_hot_in(self, fname=None):
        """ write hotstart inputs """
        with open(fname, 'w') as fout:
            fout.write(str(len(self.data['Depth'])) + '\n')
            for depth, value in zip(self.data['Depth'], self.data['MeasureValue']):
                fout.write(str(depth) + " " + str(value) + "\n")


class RUN():

    """ class for setting up model simulations """

    def __init__(self, base_run_fullpath=None, run_name=None):
        from pathlib import Path
        if base_run_fullpath is not None and run_name is not None:
            self.base_run_fullpath = os.path.abspath(base_run_fullpath)
            self.base_run_path = str(Path(self.base_run_fullpath).parent) + '/'
            self.base_run_name = os.path.basename(self.base_run_fullpath)
            print('base run path: ', self.base_run_path)
            print('base run name: ', self.base_run_name)
            self.run_name = run_name
            self.fullpath = self.base_run_path + self.run_name + '/'
            print('run name: ', self.run_name)
            # test run existence
            while (os.path.exists(self.fullpath)):
                print("Run exists: ", self.fullpath)
                input("Fix the problem then press Enter to continue...")
            else:
                print("Making new run: ", self.fullpath)
                os.makedirs(self.fullpath)
            self.file_types = ['*.gr3', '*.in', "*th*", "run.*"]
            self.nstacks = None

    def link_inputs(self):
        import glob
        for this_type in self.file_types:
            for file_from_base_run in glob.glob(self.base_run_fullpath + '/' + this_type):
                os.symlink(file_from_base_run, self.fullpath + os.path.basename(file_from_base_run))
        output_dir = '/sciclone/scr20/feiye/' + os.path.basename(self.base_run_path[:-1]) + '/' + self.run_name + '/outputs'
        if (not os.path.exists(output_dir)):
            os.system('mkdir -p ' + output_dir)
        os.symlink(output_dir, self.fullpath + '/outputs')

    def set_ops(self, sub_dir):
        import subprocess
        import shutil

        cmd = './auto.pl'
        script_dir = self.fullpath + '/' + sub_dir
        shutil.copytree(self.base_run_fullpath+'/'+sub_dir, script_dir)
        proc = subprocess.Popen(cmd, cwd=script_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = proc.communicate()

    def set_exe(self, exe_fullpath):
        exe_fname = os.path.basename(exe_fullpath)
        os.symlink(exe_fullpath, self.fullpath + '/' + exe_fname)

    def execute(self, run_cmd=None, run_arg=None):
        proc = subprocess.Popen(run_cmd, cwd=self.fullpath, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
        o, e = proc.communicate(input=run_arg)
    
    def get_para(self, var_str):
        import f90nml
        nml = f90nml.read(self.fullpath+'param.in')
        return nml['para_nml'][var_str]

    def get_nstacks(self):
        dt = self.get_para("dt")
        ihfskip = self.get_para("ihfskip")
        rnday = self.get_para("rnday")
        self.nstacks = rnday * 86400 / dt / ihfskip
    