"""
Sample usage of the source_sink module
"""

# %%
from pylib_experimental.schism_file import TimeHistory
vs1 = TimeHistory.from_file('/sciclone/schism10/feiye/STOFS3D-v7/I12w/Source_sink/relocated_source_sink2/vsource.th')
vs2 = TimeHistory.from_file('/sciclone/schism10/feiye/STOFS3D-v7/I12w/Source_sink/relocated_source_sink/vsource.th')

print(vs1 == vs2)

# %%
from pylib_experimental.schism_file import source_sink

ss = source_sink.from_ncfile(
    '/sciclone/schism10/feiye/STOFS3D-v8/I03v/03_source.nc',
)

relocate_dict = {
    131606: 173704
}

ss.reset_source_ele(relocate_dict)

ss.nc_writer('/sciclone/schism10/feiye/STOFS3D-v8/I03v/Relocated_SS/')
ss.writer('/sciclone/schism10/feiye/STOFS3D-v8/I03v/Relocated_SS/')


pass

# %%
