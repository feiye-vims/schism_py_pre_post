import os
import copy
import csv
import numpy
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
from matplotlib import pyplot
import pickle


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
