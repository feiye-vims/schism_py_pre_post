from schism_py_pre_post.Shared_modules.viz_transect import plot_snapshot


if __name__ == "__main__":
    gridfile = '/sciclone/schism10/feiye/ICOGS/RUN10g/hgrid.npz'
    bpfile = '/sciclone/schism10/feiye/ICOGS/BPfiles/missi.bp'

    # snapshot
    project = 'STOFS3D-v4'

    runids = ['RUN23k', 'RUN23k5', 'RUN23k6', 'RUN23k1', 'RUN23k4', 'RUN23k7']
    plot_snapshot(gridfile=gridfile, bpfile=bpfile, project=project, runids=runids, snapshot_idx=-1)

    # runids = ['RUN23o', 'RUN23o1']
    # plot_snapshot(gridfile=gridfile, bpfile=bpfile, project=project, runids=runids, snapshot_idx=431)
