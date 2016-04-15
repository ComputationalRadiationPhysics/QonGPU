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

def norm(a,dx):
    c = np.sum(a**2*dx)
    return 1/np.sqrt(c)*a

def R10(x):
    return x*2*np.exp(-x)


x = np.linspace(0,p[2],int(nx)-1)
dx = p[2]/nx
for i in range(0,int(len(a)/nx)-1):
    a[i:i+int(nx)-1] = norm(a[i:i+int(nx)-1],dx)
    plt.plot(x,a[i:i+int(nx)-1]**2,label="$E_{"+str(i)+"}$")
    print(np.sum(a[i:i+int(nx)-1]**2*dx))
plt.xlabel(" $x \, (a_0)$",size=20)
plt.ylabel("$|\psi_n|^2$",size=20)
l1 = norm(x*R10(x),dx)
plt.plot(x,l1**2,label = "$R_{10}$")
print(np.sum(l1**2)*dx)
plt.legend()
plt.show()
