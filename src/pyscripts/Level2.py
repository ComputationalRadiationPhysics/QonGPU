import numpy as np
import h5py
from matplotlib import pyplot as plt

f = h5py.File("static_results.h5","r")
a1 = f.get("/numres")
a2 = f.get("/params")
p = np.array(a2)
a = np.array(a1)
nx = p[0]
ne = p[1]
x = np.linspace(0,p[2],int(nx))
for i in range(0,int(len(a)/nx)-1):
    plt.plot(x,a[i:i+int(nx)-1]**2)
plt.show()
