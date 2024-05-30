import numpy as np
from matplotlib import pyplot as plt
from schism_py_pre_post.Shared_modules import extract_schism

u = np.loadtxt('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15/PostP/horizontalVelX.dat.Key_Bridge.R15')
v = np.loadtxt('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15/PostP/horizontalVelY.dat.Key_Bridge.R15')
t = u[23:, 0]; u = u[23:, 1]; v = v[23:, 1]  # remove the first 23 hours
datetime = np.array([np.datetime64('2024-03-05') + np.timedelta64(int(tt*86400), 's') for tt in t])

max_vel = np.max(np.sqrt(u**2 + v**2))
u = u / max_vel
v = v / max_vel
# add a vector to be used as the legend at the end of the time series
u = np.append(u, 0.2/max_vel)
v = np.append(v, 0)
x = np.append(np.arange(len(datetime)), len(datetime)-30)
y = np.append(np.zeros_like(t), 50)

plt.quiver(x, y, u, v, width=0.001)
plt.xticks(np.arange(0, len(datetime), 24), datetime[::24], rotation=90)
plt.yticks([])
plt.ylim([-40, 40])

# add text to the legend
plt.text(len(datetime)-30, 50, '0.2 m/s', ha='center', va='bottom')

plt.rcParams['font.size'] = 16

plt.axis('equal')
plt.tight_layout()
plt.show()

rundir = '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15/'
bpname = "Key_Bridge"
outfilenames = extract_schism(
    bpfiles=[f'/sciclone/schism10/feiye/ICOGS/BPfiles/{bpname}.bp'],
    var='horizontalVelY',
    rundir=rundir,
    i_comb=True,
    stacks=[1, 21],
    i_all_level=False,
    ver=10
)
print(f'extraction done: {outfilenames}')