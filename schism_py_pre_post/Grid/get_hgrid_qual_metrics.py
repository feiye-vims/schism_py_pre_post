import numpy as np
from pylib import read
from spp_core.Utilities.import_util import get_hgrid_reader
import sys

read_hgrid = get_hgrid_reader()
EPS = 1e-12


def area_change_ratio_from_neighbors(area, ic3):
    """
    area: (ne,)
    ic3: (ne,4) neighbor element indices, padded with -1
    Returns ACR: (ne,) = max(neigh_area)/min(neigh_area), nan if no neighbors.
    """
    neigh = ic3.copy()
    mask = neigh < 0

    neigh_area = area[neigh]
    neigh_area = neigh_area.astype(float)
    neigh_area[mask] = np.nan

    ratios = neigh_area / area[:, np.newaxis]
    ratios[ratios < 1.0] = 1.0 / np.maximum(ratios[ratios < 1.0], EPS)

    acr = np.nanmax(ratios, axis=1)
    
    return acr


def cal_skewnewss(gd):
    '''
    Calculate skewness (longest_side_length/quivalent_element_radius) of each element
    '''
    if not hasattr(gd, 'dpe'):
        gd.compute_ctr()
    if not hasattr(gd, 'distj'):
        gd.compute_side(fmt=2)
    if not hasattr(gd, 'elside'):
        gd.compute_ic3()
    if not hasattr(gd, 'area'):
        gd.compute_area()

    distj = gd.distj[gd.elside]
    distj[gd.elside == -1] = 0  # side length

    gd.skewness = distj.max(axis=1)/np.sqrt(np.maximum(gd.area/np.pi, sys.float_info.epsilon))
    return gd.skewness


def customized_histogram(
    data, bins=50, range=None, y_tick_scale=1.0,
    xlabel=None,
    ylabel='Number of Elements'
):
    """
    Customized histogram that handles nan values.
    data: 1D array
    bins: int or array-like
    range: tuple (min, max)
    Returns hist, bin_edges
    """

    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    import matplotlib as mpl
    mpl.rcParams.update({
        'font.size': 16,          # base font size
        'axes.titlesize': 16,
        'axes.labelsize': 16,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'legend.fontsize': 16,
    })

    plt.figure()
    plt.hist(data, bins=50)

    if range is not None:
        plt.xlim(range[0], range[1])

    # formatter: divide tick labels by 1000
    if y_tick_scale != 1.0:
        plt.gca().yaxis.set_major_formatter(
            FuncFormatter(lambda x, pos: f'{x/y_tick_scale:.1f}')
        )
    
    plt.xlabel(xlabel)
    plt.ylabel('Number of Elements (×{:.0f})'.format(y_tick_scale))

    # plt.title('Histogram of Element Aspect Ratios')
    plt.grid()
    plt.tight_layout()
    plt.show()


def main():
    g = read_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/R29i1/hgrid.gr3')
    g.x, g.y = g.proj(prj0='EPSG:4326', prj1='esri:102008')  # project to conus lambert
    g.compute_all()

    # aspect ratio
    distj_padded = np.r_[g.distj, np.nan]  # elside may have -1 index for none-existing sides
    element_dl = distj_padded[g.elside]
    aspect_ratios = -np.log(abs(0.5 - np.nanmax(element_dl, axis=1) / (np.nansum(element_dl, axis=1) + 1e-12)))

    # skewness
    skewness = cal_skewnewss(g)

    # element area change ratio with neighbors
    acr = area_change_ratio_from_neighbors(g.area, g.ic3)

    # plot histogram of aspect ratios
    customized_histogram(
        aspect_ratios, bins=50, range=(0, 8), y_tick_scale=1000.0,
        xlabel=r'Aspect Ratio $\left[-\log\left(0.5 - \frac{\max(\ell)}{\sum \ell}\right)\right]$',
    )

    # plot histogram of skewness
    customized_histogram(
        np.log10(skewness), bins=50, range=(0.2, 1.6), y_tick_scale=1000.0,
        xlabel='$log_{10}$ (Skewness)',
    )

    customized_histogram(
        np.log10(acr), bins=50, range=(0, 2.3), y_tick_scale=1000.0,
        xlabel='$log_{10}$ (Area Change Ratio with Neighbors)',
    )

    print('Done')

if __name__ == '__main__':
    main()