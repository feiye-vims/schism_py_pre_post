from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
import numpy as np
import os
from scipy import interpolate
import pandas as pd
import copy
from pylib import zdata, WriteNC, schism_grid


def BinA(A, B):
    import numpy as np
    x = A
    y = B

    index = np.argsort(x)
    sorted_x = x[index]
    sorted_index = np.searchsorted(sorted_x, y)

    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y

    result = np.ma.array(yindex, mask=mask)

    return [np.logical_not(np.ma.getmask(result)), np.ma.getmask(result)]


class SourceSinkIn():
    """ class for *.prop or other similar formats"""
    def __init__(self, filename=None, number_of_groups=2, ele_groups=[]):
        self.n_group = number_of_groups  # 0: source; 1: sink
        if filename is not None:
            """Read a usgs data file and initialize from its content """
            self.source_file = filename
            self.np_group = []
            self.ip_group = []

            with open(self.source_file) as fin:
                for k in range(0, self.n_group):
                    self.np_group.append(int(fin.readline().split()[0]))
                    print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                    self.ip_group.append(np.empty((self.np_group[k]), dtype=int))
                    for i in range(0, self.np_group[k]):
                        self.ip_group[k][i] = int(fin.readline())
                    fin.readline()  # empty line between groups
                    if self.np_group[k] > 0:
                        print("p first: " + str(self.ip_group[k][0]))
                        print("p last: " + str(self.ip_group[k][-1]))
        else:
            self.np_group = [len(x) for x in ele_groups]
            self.ip_group = [np.array(x) for x in ele_groups]

    def writer(self, filename=None):
        if filename is None:
            filename = self.source_file

        with open(filename, 'w') as fout:
            for k in range(0, self.n_group):
                print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                fout.write(f"{self.np_group[k]}\n")
                for i in range(0, self.np_group[k]):
                    fout.write(f"{self.ip_group[k][i]}\n")
                fout.write("\n")  # empty line

    def add_points(self, added_points=np.array([]), group=None):
        self.np_group[group] += len(added_points)
        self.ip_group[group] = np.r_[self.ip_group[group], added_points]

    def get_max(self):
        """get max prop value"""

    def get_ele(self, pval):
        """return ele id with specified prop value"""


class Prop(SourceSinkIn):  # temporary, for be backward compatible with some older scripts
    def __init__(self, filename, number_of_groups):
        super().__init__(self, filename, number_of_groups)


class source_sink():
    """class for handling all source/sink inputs"""
    def __init__(self, source_dir=None, start_time_str='2000-01-01 00:00:00', timedeltas=[0.0, 86400.0*365],
                 source_eles=[], sink_eles=[]):

        dummy_source_sink = TimeHistory(
            file_name=None, start_time_str='2000-01-01 00:00:00',
            data_array=np.c_[np.array(timedeltas), np.array(timedeltas)*0.0], columns=['datetime', '1']
        )

        if source_dir is None:
            nt = len(timedeltas)
            nsources = len(source_eles)
            nsinks = len(sink_eles)

            if nsources > 0:
                self.vsource = TimeHistory(file_name=None, start_time_str=start_time_str,
                                           data_array=np.c_[np.array(timedeltas), np.zeros([nt, nsources])],
                                           columns=['datetime']+source_eles)
                self.msource = TimeHistory(file_name=None, start_time_str=start_time_str,
                                           data_array=np.c_[np.array(timedeltas), -9999*np.ones([nt, nsources]), np.zeros([nt, nsources])],
                                           columns=['datetime']+source_eles+source_eles)
            else:
                source_eles = [1]  # dummy with 0 vsource
                self.vsource = copy.deepcopy(dummy_source_sink)
                self.msource = TimeHistory(file_name=None, start_time_str=start_time_str,
                                           data_array=np.c_[np.array(timedeltas), -9999*np.ones([nt, 1]), np.zeros([nt, 1])],
                                           columns=['datetime', '1', '1'])

            if nsinks > 0:
                self.vsink = TimeHistory(file_name=None, start_time_str=start_time_str,
                                         data_array=np.c_[np.array(timedeltas), np.zeros([nt, nsinks])],
                                         columns=['datetime']+sink_eles)
            else:
                sink_eles = [1]  # dummy with 0 vsource
                self.vsink = copy.deepcopy(dummy_source_sink)

            self.source_sink_in = SourceSinkIn(filename=None, number_of_groups=2, ele_groups=[source_eles, sink_eles])

        else:  # read from existing files
            self.source_sink_in = SourceSinkIn(f"{source_dir}/source_sink.in")
            self.source_dir = source_dir

            if self.source_sink_in.np_group[0] > 0:
                print('reading vsource\n')
                self.vsource = TimeHistory(f"{source_dir}/vsource.th", start_time_str=start_time_str, columns=['datetime']+self.source_sink_in.ip_group[0].tolist())
                print('reading msource\n')
                self.msource = TimeHistory(f"{source_dir}/msource.th", start_time_str=start_time_str, columns=['datetime']+self.source_sink_in.ip_group[0].tolist()+self.source_sink_in.ip_group[0].tolist())
                source_eles = self.source_sink_in.ip_group[0]
            else:
                self.vsource = copy.deepcopy(dummy_source_sink)
                self.msource = TimeHistory(file_name=None, start_time_str=start_time_str,
                                           data_array=np.c_[np.array(timedeltas), -9999*np.ones([nt, 1]), np.zeros([nt, 1])],
                                           columns=['datetime', '1', '1'])
                source_eles = [1]


            if self.source_sink_in.np_group[1] > 0:
                print('reading vsink\n')
                self.vsink = TimeHistory(f"{source_dir}/vsink.th", start_time_str=start_time_str, columns=['datetime']+self.source_sink_in.ip_group[1].tolist())
                self.vsink.df.columns = self.vsink.df.columns.map(str)
                sink_eles = self.source_sink_in.ip_group[1]
            else:  # dummy vsink
                self.vsink = copy.deepcopy(dummy_source_sink)
                sink_eles = [1]

            # accomodate dummy sources/sinks
            self.source_sink_in = SourceSinkIn(filename=None, number_of_groups=2, ele_groups=[source_eles, sink_eles])

        self.vsource.df.columns = self.vsource.df.columns.map(str)
        self.msource.df.columns = self.msource.df.columns.map(str)
        self.vsink.df.columns = self.vsink.df.columns.map(str)

        self.nsource = self.source_sink_in.np_group[0]
        self.nsink = self.source_sink_in.np_group[1]
        self.source_eles = source_eles
        self.sink_eles = sink_eles
        self.ntracers = 2
    
    def update_vars(self):
        self.vsource.df_propagate()
        self.msource.df_propagate()
        self.vsink.df_propagate()
        self.nsource = self.source_sink_in.np_group[0]
        self.nsink = self.source_sink_in.np_group[1]
        self.source_eles = self.source_sink_in.ip_group[0]
        self.sink_eles = self.source_sink_in.ip_group[1]
        self.ntracers = 2
        
        self.vsource.data = self.vsource.df.values[:, 1:]

    def __add__(self, other):
        A = copy.deepcopy(self)
        B = copy.deepcopy(other)
        for i, source_sink in enumerate(['vsource', 'vsink']):  # 0: source; 1: sink
            [_, new_eles_inds] = BinA(A.source_sink_in.ip_group[i],
                                      B.source_sink_in.ip_group[i])
            [existing_eles_inds, _] = BinA(B.source_sink_in.ip_group[i],
                                           A.source_sink_in.ip_group[i])
            new_eles = B.source_sink_in.ip_group[i][new_eles_inds]
            existing_eles = A.source_sink_in.ip_group[i][existing_eles_inds]

            A.source_sink_in.add_points(added_points=new_eles, group=i)

            B_ss = eval(f"B.{source_sink}")
            A_ss = eval(f"A.{source_sink}")

            f_interp = interpolate.interp1d(B_ss.time, B_ss.df.iloc[:, 1:].to_numpy().T)
            B_df_interp = pd.DataFrame(data=f_interp(A_ss.time).T, columns=B_ss.df.columns[1:])

            A_ss.df = copy.deepcopy(A_ss.df.join(B_df_interp[new_eles.astype(str)]))
            A_ss.df[existing_eles.astype('str')] += B_df_interp[existing_eles.astype('str')]

            A.msource = TimeHistory(file_name=None,
                                    data_array=np.c_[np.array([0.0, 86400*365*100]),
                                                     -9999*np.ones([2, len(A.vsource.df.columns)-1]),
                                                     np.zeros([2, len(A.vsource.df.columns)-1])],
                                    columns=['datetime']+A.vsource.df.columns[1:].tolist()+A.vsource.df.columns[1:].tolist())
            A_ss.n_time = A_ss.df.shape[0]
            A_ss.n_station = A_ss.df.shape[1] - 1

            A.update_vars()

        return A

    def writer(self, dirname=None):
        if dirname is None:
            raise Exception("dirname is required.")
        os.makedirs(dirname, exist_ok=True)

        self.vsource.writer(dirname + '/vsource.th')
        self.vsink.writer(dirname + '/vsink.th')
        self.msource.writer(dirname + '/msource.th')
        self.source_sink_in.writer(dirname + '/source_sink.in')
    
    def nc_writer(self, dirname=None):
        if dirname is None:
            raise Exception("dirname is required.")
        os.makedirs(dirname, exist_ok=True)

        # create netcdf data
        C=zdata(); C.vars=[]; C.file_format='NETCDF4'

        C.dimname=['nsources', 'nsinks', 'ntracers', 'time_msource','time_vsource','time_vsink','one']
        C.dims=[self.nsource, self.nsink, self.ntracers, self.msource.n_time, self.vsource.n_time, self.vsink.n_time, 1]

        C.vars.extend(['source_elem','vsource','msource'])
        vi=zdata(); vi.dimname=('nsources',); vi.val=self.source_sink_in.ip_group[0]; C.source_elem=vi
        vi=zdata(); vi.dimname=('time_vsource','nsources'); vi.val=self.vsource.data; C.vsource=vi
        msource_data = self.msource.data.reshape([self.msource.n_time, self.ntracers, self.nsource])
        vi=zdata(); vi.dimname=('time_msource','ntracers','nsources'); vi.val=msource_data; C.msource=vi

        if self.nsink > 0:
            C.vars.extend(['sink_elem','vsink'])
            vi=zdata(); vi.dimname=('nsinks',); vi.val=self.source_sink_in.ip_group[1]; C.sink_elem=vi
            vi=zdata(); vi.dimname=('time_vsink','nsinks',); vi.val=self.vsink.data; C.vsink=vi
        
        C.vars.extend(['time_step_vsource','time_step_msource','time_step_vsink'])
        vi=zdata(); vi.dimname=('one',); vi.val=self.vsource.delta_t; C.time_step_vsource=vi
        vi=zdata(); vi.dimname=('one',); vi.val=self.msource.delta_t; C.time_step_msource=vi
        vi=zdata(); vi.dimname=('one',); vi.val=self.vsink.delta_t; C.time_step_vsink=vi

        WriteNC(f'{dirname}/source.nc', C)


if __name__ == "__main__":
    ss = source_sink(source_dir='/sciclone/schism10/feiye/STOFS3D-v4/Inputs/IOper/')
    pass
    """testing"""
    # Sample 1: initialize source_sink files from scratch
    # added = source_sink(source_dir=None, source_eles=[23, 24, 25], sink_eles=[233, 234])
    # added.writer('/sciclone/schism10/feiye/ICOGS/RUN10a/Pump_sink_source/Test3/')

    # Sample 2: add two sets of source_sink files
    # source_sink_1 = source_sink('/sciclone/schism10/feiye/ICOGS/RUN10a/Pump_sink_source/Test1/')
    # source_sink_2 = source_sink('/sciclone/schism10/feiye/ICOGS/RUN10a/Pump_sink_source/Test2/')
    # source_sink_combined = source_sink_1 + source_sink_2
    # source_sink_combined.writer('/sciclone/schism10/feiye/ICOGS/RUN10a/Pump_sink_source/Test3/')

    # Sample 3: write source.nc in netcdf4 format
    # my_source_sink = source_sink(source_dir='/sciclone/data10/feiye/vims20/work/ChesBay/RUN200i/')
    # my_source_sink.nc_writer('./')

    # Sample 4: initialize source_sink files with uniform values
    my_hgrid = schism_grid('/sciclone/data10/feiye/vims20/work/ChesBay/RUN200i/hgrid.gr3')
    my_hgrid.compute_all()
    land_ele_ids = np.squeeze(np.argwhere(my_hgrid.dpe<5)).tolist() 
    my_sink = source_sink(source_dir=None, source_eles=[], sink_eles=land_ele_ids, timedeltas=[0, 86400*36500])
    my_source_sink = source_sink(source_dir='/sciclone/data10/feiye/vims20/work/ChesBay/RUN200i/')
    combined = my_source_sink + my_sink
    # combined.writer('./')
    combined.nc_writer('./')
    pass
