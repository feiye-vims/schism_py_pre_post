import jax.numpy as jnp
import numpy as np
import time

'''
pip install --upgrade pip

pip install --upgrade "jax[cpu]"
or
pip install --upgrade "jax[gpu]"
'''

# time step in seconds
dt = 0.005
# Frequency in Hz
f = 2.0
# Number of time steps
N = 100000000

# JAX array containing time vector consisting of 32 bit floating point precision values
t = jnp.linspace(0, N * dt, N, dtype=jnp.float32)

# Vectorized operation that passes time array into JAX sin function
jax_start_time = time.time()
signal1 = jnp.sin(2.0 * jnp.pi * f * t)
jax_end_time = time.time()

t = np.linspace(0, N * dt, N, dtype=np.float32)
np_start_time = time.time()
signal2 = np.sin(2.0 * np.pi * f * t)
np_end_time = time.time()

print(
    f"JAX CPU operation time: {jax_end_time - jax_start_time:.5f} seconds\n",
    f"Numpy operation time: {np_end_time - np_start_time:.5f} seconds",
)