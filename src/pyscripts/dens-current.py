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
    psi_t = np.zeros([nt,nx])
    print("Loading file")
    with progressbar.ProgressBar(max_value=int(nt)) as bar:
        for i  in range(0, nt, nstep):
            psi_t[i] = load_timestep(filepath, nt)
            bar.update(i)
    return psi_t
    

def conjugate(psi):
    """
    Build the complex conjugate for each timestep 
    t of psi and returns it.
    """
    psi_conj = np.conj(psi)
    return psi_conj


def get_imag_grad(psi, h):
    """
    Numerically differentiates the complex conjugated
    wave-function psi with the spatial step h. The parameter
    interval should be a tuple of indexes, from when to where to
    differentiate from.  
    """
    psi_conj = conjugate(psi)
    psi_diff = np.gradient(psi_conj, h, axis=1)
    return psi_diff

def get_prob_current(psi, psi_diff):
    """
    Calculates the probability current Im{psi d/dx psi^*}. 
    """
    curr = psi*psi_diff
    return curr.imag


def integrate_prob_current(psi, x0, x1, h):
    """
    Numerically integrate the probability current, which is
    Im{psi d/dx psi^*} over the given spatial interval. 
    """
    psi_diff = get_imag_grad(psi, h)
    curr = get_prob_current(psi, psi_diff)
    
    # Using an int32 as indexes
    n0 = np.int32(x0/h)
    n1 = np.int32(x1/h)
    res = np.zeros(psi.shape[0])
    with progressbar.ProgressBar(max_value=int(psi.shape[0])) as bar:
        for i in range(0, psi.shape[0]):
            res [i] = np.trapz(curr[i,n0:n1], dx=h)
            bar.update(i)
    
    return res

def main():

    filepath = "../../simres/sim_axel.h5"
    nx = 1e5
    nt = 1e5
    nstep = 100
    psi = load_vals(filepath, nt, nx, nstep)

if __name__ == "__main__":
    main()
