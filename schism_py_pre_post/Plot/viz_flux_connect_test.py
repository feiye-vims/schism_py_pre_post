"""
Visualize channel connectivity by comparing upstream fluxes with downstream fluxes.
Ideally, the upstream fluxes should be equal to the downstream fluxes.
For each group (selected watershed), the first transect is the downstream transect,
and the remaining transects are upstream transects, which should all flow into the first transect.
"""
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory

groups = {
    'a': ['a00', 'a01', 'a02', 'a03', 'a04', 'a05', 'a06', 'a07', 'a08'],
    # 'b': ['b00', 'b08', 'b09', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15'],
    'b': ['Downstream', '1', '2', '3', '4', '5', '6', '7', '8'],
    'c': ['c00', 'c16', 'c17', 'c18', 'c25'],
    'd': ['d00', 'd01', 'd02', 'd03', 'd04', 'd05', 'd06', 'd07', 'd08', 'd09', 'd10'],
    'e': ['e00', 'e01', 'e03'],
    'f': ['f00', 'f02'],
}
columns = ['datetime'] + [i for group, ids in groups.items() for i in ids]

flux = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v8/R09c2/outputs/flux.out',
                   start_time_str="2024-03-05 00:00:00", sec_per_time_unit=86400)
# time_idx = np.argwhere(np.array(flux.df['datetime'] >= datetime(2024, 3, 20, 0, 0, 0))).flatten()
# flux = flux.export_subset(time_idx=time_idx)
flux.df.columns = columns

# # special case for 'e' and 'd'
# downstream_transects = ['e00', 'f00']
# upstream_transects = [transect for transect in groups['e'][1:] + groups['f'][1:]]
# total_upstream = np.sum(flux.df[upstream_transects], axis=1)
# total_downstream = np.sum(flux.df[downstream_transects], axis=1)
# 
# fig = plt.figure(figsize=(12, 7))
# for transect in upstream_transects + downstream_transects:
#     plt.plot(flux.df['datetime'], flux.df[transect], label=transect)
# plt.plot(flux.df['datetime'], total_upstream, label='total upstream', color='black', linestyle='--')
# plt.plot(flux.df['datetime'], total_downstream, label='total downstream', color='red', linestyle='--')
# 
# print(f'percentage: {total_downstream.mean() / total_upstream.mean() * 100:.2f}%')
# plt.legend()
# plt.show()

plot_groups = {k: groups[k] for k in ['b'] if k in groups}

plt.rcParams.update({'font.size': 20})
fig, axs = plt.subplots(len(plot_groups) + 1, 1, figsize=(12, 7), sharex=True)
for i, (group, transects) in enumerate(plot_groups.items()):
    axs[i+1].plot(
        flux.df['datetime'], flux.df[transects[0]], label=transects[0],
        color='#5fa2d0', linestyle='-', linewidth=4
    )  # downstream

    for transect in transects[1:]:  # individual upstream transects
        axs[i+1].plot(flux.df['datetime'], flux.df[transect], label=f'Transect {transect}')
    total_flow = np.sum(flux.df[transects[1:]], axis=1)

    axs[i+1].plot(
        flux.df['datetime'], total_flow, label='Sum of upstream flow',
        color='black', linestyle='--', linewidth=4
    )  # total upstream flow

    percentage = flux.df[transects[0]].mean() / np.mean(total_flow) * 100
    title_str = f'{group} {transects[0]} = {percentage:.2f}% of total upstream flow'

    axs[i+1].set_title(title_str)
    axs[i+1].legend()
    axs[i+1].set_ylabel(r'Flow ($m^3\ s^{-1}$)')

plt.xticks(rotation=25)
plt.show()
print('Done')
