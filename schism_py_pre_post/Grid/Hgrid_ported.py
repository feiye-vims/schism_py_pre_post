from pylib import schism_grid
import os
import numpy as np
import pickle
import copy

from sqlalchemy import over


def read_schism_hgrid_cached(gd_filename, overwrite_cache=False):

    # gd_cache_fname = os.path.splitext(gd_filename)[0] + '.pkl'
    gd_cache_fname = gd_filename + '.pkl'

    if overwrite_cache or not os.path.exists(gd_cache_fname):
        gd = schism_grid(gd_filename)
        with open(gd_cache_fname, 'wb') as file:
            pickle.dump(gd, file)
    else:
        with open(gd_cache_fname, 'rb') as file:
            gd = pickle.load(file)

    return gd

