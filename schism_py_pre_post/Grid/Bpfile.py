import numpy as np


class Bpfile():

    """class for bpfiles"""

    def __init__(self, filename, cols=4):
        """Initialization """
        self.nodes = []
        self.n_nodes = 0
        self.info = []

        self.reader(filename, cols)

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
        self.nodes = np.array(self.nodes, dtype=float)

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
            for i, st in enumerate(self.info):
                self.st_id[i] = st[0][1:]
            col_name = self.info

        self.df = pd.DataFrame(data=self.nodes.T, index=row_name, columns=col_name)
        return self.df
