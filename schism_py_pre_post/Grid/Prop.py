import numpy 


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



