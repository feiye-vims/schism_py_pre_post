
import numpy as np
import matplotlib.pyplot as plt
from pylib_experimental.schism_file import TimeHistory, SourceSinkIn


def main():
    ssin = SourceSinkIn.from_file('/sciclone/schism10/feiye/STOFS3D-v8/I15n_v7/Source_sink/original_source_sink/source_sink.in')
    th1 = TimeHistory.from_file('/sciclone/schism10/feiye/STOFS3D-v8/I15n_v7/Source_sink/USGS_adjusted_sources/adjusted_vsource.th', start_time_str='2017-12-01 00:00:00')
    th2 = TimeHistory.from_file('/sciclone/schism10/feiye/STOFS3D-v7.3/I18/Source_sink/USGS_adjusted_sources/adjusted_vsource.th', start_time_str='2017-12-01 00:00:00')
    # assert th1 == th2

    diff = th1.data - th2.data
    mean_diff = np.mean(np.abs(diff), axis=0)
    large_diff_idx = np.argwhere(mean_diff > 0.1).flatten()
    fig, ax = plt.subplots(ncols=2, nrows=len(large_diff_idx)//2+1, figsize=(12, 4))
    for i, idx in enumerate(large_diff_idx):
        ele = ssin.ele_groups[0][idx]
        ax[i//2, i%2].set_title(f'Element {ele} Index {idx} Mean Abs Diff: {mean_diff[idx]:.6f}')
        ax[i//2, i%2].plot(th1.datetime, th1.data[:, idx], label='th1')
        ax[i//2, i%2].plot(th2.datetime, th2.data[:, idx], '--r', label='th2')
        ax[i//2, i%2].legend()
    plt.show()

    plt.plot(mean_diff)
    plt.show()
    print("Max diff:", max(diff.flatten()))


if __name__ == "__main__":
    main()