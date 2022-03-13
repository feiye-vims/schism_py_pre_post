import numpy as np


class Bpfile():

    """class for bpfiles"""

    def __init__(self, filename=None, cols=4, nodes=None):
        """Initialization """
        self.nodes = []
        self.n_nodes = 0
        self.info = []

        if filename is not None:
            self.reader(filename, cols)
        else:
            nodes = np.array(nodes)
            self.n_nodes = len(nodes[:, 0])
            if len(nodes[:, 1] < 3):
                nodes = np.c_[nodes, np.zeros([self.n_nodes, 1])]
            self.nodes = nodes
            for i in range(self.n_nodes):
                self.info.append(f"Station_{i+1}")

        print(f"{self.n_nodes}")
        print(f"{self.nodes[0][:]}")

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

    def writer(self, out_file_name):
        """ nodes should be shaped as [:n_nodes][:2] """
        with open(out_file_name, 'w') as fout:
            fout.write("\n")
            fout.write(f"{self.n_nodes}\n")
            for i in range(self.n_nodes):
                fout.write(f"{i+1} {self.nodes[i][0]} {self.nodes[i][1]} {self.nodes[i][2]} !{self.info[i]}\n")

    def make_dataframe(self, row_name=['lon', 'lat', 'z'], col_name=None):
        import pandas as pd
        if col_name is None:
            for i, st in enumerate(self.info):
                self.st_id[i] = st[0][1:]
            col_name = self.info

        self.df = pd.DataFrame(data=self.nodes.T, index=row_name, columns=col_name)
        return self.df
