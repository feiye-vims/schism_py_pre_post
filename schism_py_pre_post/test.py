# %%
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
%matplotlib inline
set_matplotlib_formats('svg', 'png', 'pdf')

plt.plot([1, 2, 3], [34, 343, 3423])
plt.savefig("test.svg", format="svg")
# %%
