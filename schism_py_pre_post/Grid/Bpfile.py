import numpy as np


class Bpfile():

    """class for bpfiles"""

    def __init__(self, filename=None, cols=4, xyz_array=None):
        """Initialization """
        self.nodes = []
        self.n_nodes = 0
        self.info = []
        self.ncol = 0

        if xyz_array is not None:
            self.nodes = xyz_array
            self.n_nodes = xyz_array.shape[0]
            self.ncol = xyz_array.shape[1] + 1
        else:
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
            self.st_id = []
            if self.info != []:
                for i, st in enumerate(self.info):
                    if self.ncol < 4:
                        self.st_id.append(f'{i}')
                    else:
                        if len(st)>0:  # station name
                            self.st_id.append(st[0][1:])
                        else:
                            self.st_id.append(f'Station_{i+1}')
            else:
                self.st_id = [f'Station_{i+1}' for i in range(self.n_nodes)]
            col_name = self.st_id

        self.df = pd.DataFrame(data=self.nodes.T, index=row_name, columns=col_name)
        return self.df

