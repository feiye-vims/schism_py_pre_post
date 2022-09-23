from pylib import loadz, npz_data, save_npz
from pandas import read_csv, to_numeric
import pickle
import os
import glob

def pickle_save(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def pickle_load(filename):
    """ Deserialize a file of pickled objects. """
    with open(filename, "rb") as f:
        data = pickle.load(f)
        return data

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
        self.id = ''
        self.lon = -999
        self.lat = -999
        self.df = None

        if fname is not None:
            self.id = os.path.basename(fname).split('_')[1]
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
