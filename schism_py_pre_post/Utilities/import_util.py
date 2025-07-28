"""
This module contains utility functions for importing modules and packages.
"""


import socket


def get_hgrid_reader():
    """
    Get the appropriate schism reader based on the hostname.
    """

    if 'gulf' in socket.gethostname():
        from pylib_experimental.schism_file import cread_schism_hgrid as read_hgrid
        print('Using c++ function to accelerate hgrid reading')
    else:
        from pylib import schism_grid as read_hgrid
        print('Using python function from pylib to read hgrid')
    
    return read_hgrid