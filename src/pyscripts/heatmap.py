
from __future__  import division, print_function

import matplotlib as mpl




import matplotlib.pyplot as plt
import numpy as np
import h5py

print(mpl.__version__)


file_loc = "../../simres/heat_res.h5"
file = h5py.File(file_loc)

str_0 = "/dset"

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


plt.rc('font', family='serif')
plt.rc('xtick', labelsize=30)
plt.rc('ytick', labelsize=30)


t = np.linspace(0, 5000, 1000)
x = np.linspace(-30,30,100000)

#T,X = np.meshgrid(t,x,sparse=False,indexing='xy')

pic = np.zeros((100000,1000))
#prob = np.zeros(1000) 
dimag = 0
dreal = 0
psi = 0
print("Working")
for i in range(0, 100000, 100):
    
    ind = i+1
    
    r_str = str_0 + str(ind) + "real"
    i_str =str_0 + str(ind) + "img"
    
    if r_str in file and i_str in file:
        
        dreal = np.array(file[r_str])
        dimag = np.array(file[i_str])
        psi = dreal**2 + dimag**2
        dreal = 0
        dimag = 0
        
    else:
        print("Not found!")
        file.close()
        exit()
    #prob[int(i/100)] = np.trapz(y=psi,x=x)
    #print(np.trapz(y=psi,x=x))
    pic[:,int(i/100)] = psi


print("Plotting!")


g_t = lambda t: np.exp(-a*(t-t0)**2)
g_s = lambda x: np.exp(-b*(x)**2)

a = np.log(2)/(1500-2500)**2
print("a = "+str(a))
b = np.log(2)/(10)**2
x = np.linspace(-30, 30, 1e3)

t0 = 2500
x0 = 0

plt.imshow(pic,cmap="magma_r",extent=[t[0],t[-1],x[0],x[-1]],aspect='auto',interpolation='nearest')
plt.plot(t,10*g_t(t)-15, color="white", linewidth=3)
plt.text(510,-12,r"$g_t(t) $",fontsize=30,color="white")
plt.xlabel(r"$t  \ (a.u.)$",size=40)
plt.ylabel(r"$x \ (a.u.)$",size=40)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
print("Plotting finished!")
c = plt.colorbar()
c.set_label(r"Probability density  $|\psi_n|^2$",size=30)
print("Colorbar done!")
plt.show()
print("Showing done!")
file.close()
print("Finised!")
