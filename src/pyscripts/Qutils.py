import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import progressbar
import sys


def load_timestep(filepath,nt,t0=0,):
    """
    Load a timestep from a result file
    and returns it as a complex numpy array.
    """

    file = h5py.File(filepath)
    rl = "/dset"+str(nt)+"real"
    im = "/dset"+str(nt)+"img"
    if (rl in file) and (im in file):
        imag = np.array(file[rl])
        real = np.array(file[im])
        res = real + 1j * imag
        file.close()
        return res
    else:
        print("Error timestep "+str(nt)+" not found!")
        file.close()
        sys.exit("Data does not exist or is corrupted")
        return -1


def load_vals(filepath, nt, nx, nstep):
    """
    Loads the whole time dependent wave-function in memory.
    """
    psi_t = np.zeros([np.int32(nt/nstep),nx], dtype=complex)
    print("Loading file")
    with progressbar.ProgressBar(max_value=int(nt)) as bar:
        for i  in range(0, nt, nstep):
            psi_t[int(i/nstep)] = load_timestep(filepath, i)
            bar.update(i)
    return psi_t


def get_imag_grad(psi, h):
    """
    Numerically differentiates the complex conjugated
    wave-function psi with the spatial step h. The parameter
    interval should be a tuple of indexes, from when to where to
    differentiate from.  
    """
    print("Calculating conjugate")
    psi_conj = np.zeros(psi.shape, dtype=complex)
    with progressbar.ProgressBar(max_value=int(psi.shape[0])) as bar:
        for i in range(0, psi.shape[0]):
            psi_conj[i] = np.conjugate(psi[i,:])
            bar.update(i)
    print("Calculating gradient...")
    psi_diff = np.gradient(psi_conj, h, axis=1)
    print("Finished gradient!")
    return psi_diff

def get_prob_current(psi, psi_diff):
    """
    Calculates the probability current Im{psi d/dx psi^*}. 
    """
    print("Calculating probability current")
    curr = psi*psi_diff
    return -curr.imag


def integrate_prob_current(psi, n0, n1, h):
    """
    Numerically integrate the probability current, which is
    Im{psi d/dx psi^*} over the given spatial interval. 
    """
    psi_diff = get_imag_grad(psi, h)
    curr = get_prob_current(psi, psi_diff)
    
    res = np.zeros(psi.shape[0])
    with progressbar.ProgressBar(max_value=int(psi.shape[0])) as bar:
        for i in range(0, psi.shape[0]):
            res [i] = np.trapz(curr[i,n0:n1], dx=h)
            bar.update(i)
    print("Finished calculating the integrated prob. current!")
    return res

def integrate_prob(psi, n0, n1, h):
    """
    Integrate the |psi|^2
    """
    res = np.zeros(psi.shape[0])
    with progressbar.ProgressBar(max_value=int(psi.shape[0])) as bar:
        for i in range(0, psi.shape[0]):
            res[i] = np.trapz(np.abs(psi[i,n0:n1])**2,dx=h)
            bar.update(i)
    return res


def get_prob_dens(psi, nt, nstep):
    """
    Calculates the probability density for
    a given wavefunction psi recorded in nt
    timesteps.
    """
    res = np.zeros(psi.shape)
    for i in range(0, psi.shape[0]):
        res[i] = psi[i,:].imag**2+psi[i,:].real**2
    return res
