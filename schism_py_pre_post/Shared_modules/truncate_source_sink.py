from schism_py_pre_post.Grid.SourceSinkIn import source_sink


ss = source_sink(source_dir='/sciclone/schism10/feiye/ICOGS/RUN22/', start_time_str='2021-08-10 00:00:00')
tidx = ss.vsource.df['datetime'] > '2021-08-09 23:59:59'
ss.vsource.export_subset(station_idx=[], time_idx=tidx, i_reset_time=1, subset_filename='/sciclone/schism10/feiye/ICOGS/RUN22/Truncate_ss/vsource.th')
ss.vsink.export_subset(station_idx=[], time_idx=tidx, i_reset_time=1, subset_filename='/sciclone/schism10/feiye/ICOGS/RUN22/Truncate_ss/vsink.th')

ss_trun = source_sink(source_dir='/sciclone/schism10/feiye/ICOGS/RUN22/Truncate_ss/', start_time_str='2005-08-10 00:00:00')
pass
