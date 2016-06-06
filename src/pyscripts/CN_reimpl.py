import numpy as np
import h5py
import scipy.linalg as lag
from matplotlib import pyplot as plt



def pot(x):
    return -1/np.sqrt(x**2 + 1)

file = h5py.File("../../build/copytest.h5")



x = np.array(file["/dset0real"])

file.close()

file = h5py.File("../../build/res.h5")
params = np.array(file["/params"])
file.close()

y = np.zeros(x.size)

z = x + 1j*y


print(params)


x,h = np.linspace( params[1], params[0], params[4], retstep=True)
t,tau = np.linspace( params[2], params[3], params[5], retstep=True)
k1 = 1j*-tau/(4*h**2)*np.ones(x.size-1)

H1 = np.diag(1 + 1j * tau/2*( z/h**2 + pot(x)))
H2 = np.diag(-k1,k=-1)
H3 = np.diag(-k1, k=1)
HL = H1 + H2 + H3

