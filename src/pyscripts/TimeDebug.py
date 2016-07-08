# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
import h5py




def  V(x, t):

    return -2/(x**2+1)

def mptridiag(psi, psi2, tau, h, x):

    d1 = tau/2
    d2 = -1/(2*h**2)

    psi[0] = psi2[0] - 1j * d1 * (  d2 * ( psi2[1] - 2*psi2[0]) +  V(0,0))
    psi[-1] = psi2[-1] - 1j * d1 *( d2 * (psi2[-2] - 2 * psi2[-1]) + V(x[-1],0))
    print(psi[0])
    print(psi[-1])
    for i in range(1, (psi2.size-2)):

        psi[i] = psi2[i] - 1j * d1 * ( d2 * ( psi2[i+1] + psi2[i-1] - 2 * psi2[i]) + V(x[i],0))



file = h5py.File("../../build/res.h5")

dsetimg = np.array(file["/dset0img"])
dsetrl = np.array(file["/dset0real"])
params = np.array(file["/params"])
print(params)

xmax = params[0]
xmin = params[1]
tmax = params[2]
tmin = params[3]
nx   = params[4]
nt   = params[5]
en   = params[7]

x,h = np.linspace(xmin, xmax, nx, retstep=True)
tau = (tmax)/nt
print(1/(-2*h**2),tau)

fig = plt.figure()
ax = fig.add_subplot(111)

psi1 = dsetrl + 1j * dsetimg

psi2 = psi1

mptridiag(psi1, psi2, tau, h, x)

ax.plot(psi1.imag)
print(psi1[0:100])
plt.show()

file.close()