#/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import h5py
import numpy as np
import functools
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
#from sympy.physics.hydrogen import R_nl
#from sympy import var
import gc


f = h5py.File("static_results.h5",'r')
dset1 = f.get("/numres")
dset2 = f.get("/params1d")
dset3=f.get("/levels")
p = np.array(dset2)
a = np.array(dset1)
ind=np.array(dset3)
print(len(ind))
f.close()
dset1=0;
dset2=0;
dset3=0
nx=int(p[4])
ne = int(p[6])

x=np.linspace(0,p[0],nx-1)
dx=p[0]/nx
fig=plt.figure()
ax=fig.add_subplot(111)
# a[int(ind[i])*nx:int(ind[i])*nx+nx] = normalize(a[int(ind[i])*nx:int(ind[i])*nx+nx],p[0],nx)

def norm(a,dx):
    c = np.sum(a**2*dx)
    return 1/np.sqrt(c)*a

def R10(x):
    br = 5.29e-10
    return x*2*np.exp(-x)

R20= lambda r:r*np.sqrt(2)*(-r + 2)*np.sqrt(1**3)*np.exp(-r/2)/4
R30= lambda r:r*2*np.sqrt(3)*(2*r**2/9 - 2*1*r + 3)*np.sqrt(1**3)*np.exp(-1*r/3)/27

for i in range(0,len(ind)): 
    a[int(ind[i])*nx+1:int(ind[i])*nx+nx] = norm(a[int(ind[i])*nx+1:int(ind[i])*nx+nx],dx)
    ax.plot(x,a[int(ind[i])*nx+1:int(ind[i])*nx+nx]**2,label="$E_{"+str(i)+"}$")
#    ax.plot(x,a[i*nx+1:(i+1)*nx])

#ax.plot(x,R20(x)**2,label="$2s$-Orbital")
#ax.plot(x,R30(x)**2,label="$3s$-Orbital")
ax.set_xlabel("$x$ $(a_0)$",size=20)
ax.set_ylabel("$|\Psi_n|^2$",size=20)
#plt.legend()
i=0
plt.legend()
plt.show()
x=0
a=0
p=0
