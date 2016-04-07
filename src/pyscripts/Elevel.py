#/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import h5py
import numpy as np
import functools
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from sympy.physics.hydrogen import R_nl
from sympy import var
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
for i in range(0,len(ind)): 
    ax.plot(x,a[int(ind[i])*nx+1:int(ind[i])*nx+nx],label="$E_{"+str(i)+"}$")
#    ax.plot(x,a[i*nx+1:(i+1)*nx])
ax.set_xlabel("$x$ $(a_0)$",size=20)
ax.set_ylabel("$|\Psi_n|^2$ (not normed)",size=20)
#plt.legend()
i=0
plt.legend()
plt.show()
x=0
a=0
p=0
