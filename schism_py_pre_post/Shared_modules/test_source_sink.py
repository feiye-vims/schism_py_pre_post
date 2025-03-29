"""
Sample usage of the source_sink module
"""

from pylib_experimental.schism_file import TimeHistory
from pylib_experimental.schism_file import source_sink


def diff_vsource():
    """
    Check if two vsource files are different
    """
    vs1 = TimeHistory.from_file('/sciclone/schism10/feiye/STOFS3D-v7/I14x/Source_sink/vsink.th')
    vs2 = TimeHistory.from_file('/sciclone/schism10/feiye/STOFS3D-v7/I14x/00x_ss/vsink.th')

    assert vs1 == vs2, 'Two vsource files are different'


def manual_relocate():
    """
    Manually relocate the source elements
    """

    ss = source_sink.from_files(
        '/sciclone/schism10/feiye/STOFS3D-v8/I09/Source_sink/USGS_adjusted_sources/',
    )

    relocate_dict = {
        63209: 112941
    }

    ss.reset_source_ele(relocate_dict)

    ss.nc_writer('/sciclone/schism10/feiye/STOFS3D-v8/I09/Relocated_SS/')
    ss.writer('/sciclone/schism10/feiye/STOFS3D-v8/I09/Relocated_SS/')


if __name__ == '__main__':
    diff_vsource()
    manual_relocate()
    print('Done')
