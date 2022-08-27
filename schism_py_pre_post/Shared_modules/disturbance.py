# %%
from math import dist
from pylib import schism_grid
import copy
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from schism_py_pre_post.Plot.manual_cmap import make_colormap
from IPython.display import set_matplotlib_formats

set_matplotlib_formats('png', quality=200)

disturb_thres = 0.0
dominant_ratio = 0.8

def cal_disturb(hg, maxelev, land_node_idx, h0=1e-5, disturb_thres=0.0):
    # in ocean
    disturb = copy.deepcopy(maxelev)
    # onland
    disturb.dp[land_node_idx] = maxelev.dp[land_node_idx] + hg.dp[land_node_idx]
    # abs ?
    disturb.dp[:] = np.abs(disturb.dp[:])
    # mask dry points
    dry = maxelev.dp + hg.dp < h0
    disturb.dp[dry] = 0.0
    # mask small disturbance
    disturb.dp[disturb.dp < disturb_thres] = 0.0

    return disturb

# gd = schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10d1/PostP/elev_max_stack121-121.pkl')
# gd[0].save('/sciclone/schism10/feiye/ICOGS/RUN10d1/PostP/elev_max_stack121-121.pkl')

hg = schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10d/PostP/hgrid.pkl')  # hg.save('/sciclone/schism10/feiye/ICOGS/RUN10d/PostP/hgrid.pkl')
hg_incity = schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10e2/PostP/in_city.pkl')
incity = hg_incity.dp == 1
onland = hg.dp < 0
dry_land = (onland | incity)  # in watershed or city

maxelevs = [
    schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10b/PostP/elev_max_stack121-121.pkl'),  # ocean
    schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10d1/PostP/elev_max_stack121-121.pkl'),  # river
    schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10e2/PostP/elev_max_stack121-121.pkl'),  # precipitation
]
quiescent_disturb = cal_disturb(
    hg, schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10g/PostP/elev_max_stack121-121.pkl'),  # nothing
    dry_land, disturb_thres=disturb_thres
)
baseline_disturb = cal_disturb(
    hg, schism_grid('/sciclone/schism10/feiye/ICOGS/RUN10a/PostP/elev_max_stack121-121.pkl'),  # nothing
    dry_land, disturb_thres=disturb_thres
)
baseline_disturb.dp[:] -= quiescent_disturb.dp[:]

disturbs = []
total_disturb = copy.deepcopy(maxelevs[0])
total_disturb.dp[:] = 0
for i, maxelev in enumerate(maxelevs):
    disturb = cal_disturb(hg, maxelev, dry_land, disturb_thres=disturb_thres)
    # remove the quiescent elev
    disturb.dp = disturb.dp - quiescent_disturb.dp

    disturbs.append(disturb)
    total_disturb.dp += disturb.dp

trivial_disturb = total_disturb.dp < disturb_thres
total_disturb.dp[trivial_disturb] = np.nan

idx = total_disturb.dp > -99999
disturb_ratios = copy.deepcopy(maxelevs)
for disturb_ratio, disturb in zip(disturb_ratios, disturbs):
    disturb_ratio.dp[idx] = disturb.dp[idx] / (np.maximum(total_disturb.dp[idx], 1e-20))
    disturb_ratio.dp[~idx] = 0.0

dominance = copy.deepcopy(maxelevs[0])
dominance.dp[:] = 0  # compound
for i, disturb_ratio in enumerate(disturb_ratios):
    idx = disturb_ratio.dp >= dominant_ratio  # > dominant_ratio is considered as dominant
    dominance.dp[idx] = i + 1
dominance.dp[trivial_disturb] = np.nan

# %%
# plot disturb_ratio
plt.rcParams.update({'font.size': 20})
fig, ax = plt.subplots(3, 1, figsize=(14, 28))
# ax = ax.ravel()
# for i, disturb in enumerate(disturbs):
for i, gd in enumerate([quiescent_disturb, baseline_disturb, total_disturb]):
    plt.subplot(len(disturb_ratios), 1, i + 1)
    gd.plot_grid(fmt=1, clim=[0, 1], levels=100, ticks=np.linspace(0.0, 1, 9), cmap='jet')
    plt.xlim([-92.5, -89])
    plt.ylim([29, 31])
    plt.gca().set_aspect('equal', 'box')
plt.show()
plt.savefig('disturb.png', dpi=600)

pass
# %%
# plot disturb_ratio
plt.rcParams.update({'font.size': 20})
fig, ax = plt.subplots(len(disturb_ratios), 1, figsize=(14, 28))
# ax = ax.ravel()
# for i, disturb in enumerate(disturbs):
for i, gd in enumerate(disturb_ratios):
    plt.subplot(len(disturb_ratios), 1, i + 1)
    gd.plot_grid(fmt=1, clim=[0, 1], levels=100, ticks=np.linspace(0.0, 1, 9), cmap='jet')
    plt.xlim([-92.5, -89])
    plt.ylim([29, 31])
    plt.gca().set_aspect('equal', 'box')
plt.show()
plt.savefig('disturb.png', dpi=600)

pass

# %%
# plot dominance
plt.rcParams.update({'font.size': 30})
fig, ax = plt.subplots(1, 1, figsize=(20, 12))
c = mcolors.ColorConverter().to_rgb
gbry = make_colormap([
    c('grey'), c('grey'), 0.1,
    c('grey'), c('blue'), 0.3,
    c('blue'), c('red'), 0.5,
    c('red'), c('yellow'), 0.7,
    c('yellow'), c('yellow'), 0.9,
    c('yellow'),
])
dominance.plot_grid(fmt=1, clim=[0, 4], levels=6, ticks=np.arange(0, 4), cmap=gbry)
plt.xlim([-92.5, -89])
plt.ylim([29, 31])
plt.gca().set_aspect('equal', 'box')
plt.show()
plt.savefig('dominance.png', dpi=400)