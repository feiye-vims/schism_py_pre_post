import numba
import numpy as np
from numba.typed import List
import time

@numba.jit(nopython=True)
def calculate_numba(T, gamma):
    Nt = len(T)
    Nx = len(T[0])
    for k in range(0, Nt - 1, 1):
        for i in range(1, Nx - 1, 1):
            for j in range(1, Nx - 1, 1):
                T[k + 1, i, j] = gamma * (T[k, i + 1, j] + T[k, i - 1, j] + T[k, i, j + 1] + T[k, i, j - 1])
    return T

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    L = 50
    Nt = 1000
    Nx = 50
    alpha = 2.0
    dx = L/Nx
    dt = dx**2/4.0/alpha
    gamma = alpha*dt/dx/dx
    T_top = 100.0
    T_left = 0.0
    T_right = 0.0
    T_bottom = 0.0
    T_initial = 0.0

    # Initialize Numpy array T to store temperature values
    T = np.full((Nt,Nx,Nx),T_initial,dtype=float)
    T[:,:,:1] = T_left
    T[:,:,Nx-1] = T_right
    T[:,:1,:] = T_bottom
    T[:,Nx-1:, :] = T_top

    start_time = time.time()
    T_numba = calculate_numba(T, gamma)
    print(f"Numba JIT (first run): = {time.time() - start_time:.3f} seconds")
    start_time = time.time()
    T_numba = calculate_numba(T, gamma)
    print(f"Numba JIT (second run): = {time.time() - start_time:.3f} seconds")

    fig, ax = plt.subplots(ncols=3,figsize=(8,3),sharey='row')
    k = 999
    data = {'Python loops':T_numba[k], 'Numpy vectorized':T_numba[k], 'Numba JIT':T_numba[k]}
    i = 0
    for key, value in data.items():
        pcm = ax[i].pcolormesh(value, cmap=plt.cm.viridis, vmin=0, vmax=100)
        ax[i].set_xlabel('x-position')
        ax[i].set_aspect('equal')
        ax[i].annotate(key, xy=(1,1), c='white', fontsize=12)
        fig.colorbar(pcm,ax=ax[i],shrink=0.75)
        i+=1
    ax[0].set_ylabel('y-position')
    #fig.colorbar(pcm)
    plt.tight_layout()
    plt.show()

