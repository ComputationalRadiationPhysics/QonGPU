import numpy as np
import h5py
from matplotlib import pyplot as plt

f = h5py.File("res.h5","r")
a1 = f.get("/numres")
a2 = f.get("/params")
p = np.array(a2)
a = np.array(a1)
nx = p[0]
ne = p[1]
for i in range(0,int(len(a))/nx):
    plt.plot(a[i:i+nx])
plt.show()
