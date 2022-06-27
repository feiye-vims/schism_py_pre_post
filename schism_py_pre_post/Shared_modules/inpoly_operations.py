import numpy as np


def set_val(vals_in, idx, val_spec=None):
    if hasattr(vals_in, "__len__"):
        vals_in[:] = val_spec
    else:
        vals_in = val_spec
    return vals_in


def set_min(vals_in, idx, val_spec=None):
    vals_in = np.minimum(vals_in, val_spec)
    return vals_in


def set_max(vals_in, idx, val_spec=None):
    vals_in = np.maximum(vals_in, val_spec)
    return vals_in


def add_val(vals_in, idx, val_spec=None, val_thres=None):
    if val_thres is None:
        vals_in = (vals_in + val_spec)
    else:
        vals_in = np.minimum((vals_in + val_spec), val_thres)
    return vals_in


def subtract_val(vals_in, idx, val_spec=None, val_thres=None):
    if val_thres is None:
        vals_in = (vals_in - val_spec)
    else:
        vals_in = np.maximum((vals_in - val_spec), val_thres)
    return vals_in
