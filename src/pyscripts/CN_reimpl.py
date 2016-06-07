import numpy as np
import h5py
import scipy.linalg as lag
from matplotlib import pyplot as plt



def pot(x):
    return x**2/2

# file = h5py.File("../../build/copytest.h5")



# x1 = np.array(file["/dset0real"])

# file.close()

# file = h5py.File("../../build/res.h5")
# params = np.array(file["/params"])
x, h = np.linspace(-6, 6, 1e3, retstep=True)

# file.close()

x1 = 2*x*1/np.pi**(0.25)*1/np.sqrt(2)*np.exp(-x**2/2)
y = np.zeros(x1.size)

z = x1 + 1j * y


# print(params)


t,tau = np.linspace( 0, 10, 1e4, retstep=True)
k1 = 1j*-tau/(4*h**2)*np.ones(x1.size - 1)

H1 = np.diag(1 + 1j * tau / 2 * (1 / h**2 + pot(x)))
H2 = np.diag(k1,k=-1)
H3 = np.diag(k1, k=1)
HL = H1 + H2 + H3


H1 = np.diag(1 - 1j * tau / 2 * (1 / h**2 + pot(x)))
H2 = np.diag(-k1,k=-1)
H3 = np.diag(-k1, k=1)
HR = H1 + H2 + H3


numt = 30


def U(En, nt, tau):

    return np.exp(-1j * En * numt * tau)

psi_t = np.dot(HR, z)
fig = plt.figure()
ax = fig.add_subplot(111)
p = ax.plot(x,psi_t.imag)
an = z
for i in range(0, numt):

    print(U(1.5, i+1, tau))
    psi_t = lag.solve(HL, psi_t)
    plt.setp(p[0],ydata=psi_t.imag)
    ax.relim()
    ax.autoscale_view(True, True, True)
    plt.pause(0.05)
    psi_t = np.dot(HR, psi_t)


an= U(1.5, numt, tau) * z
plt.plot(x,an.imag)
plt.draw()
plt.show()