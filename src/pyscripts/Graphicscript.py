import numpy as np
from matplotlib import pyplot as plt
import h5py
import functools

i=0
def drawsol(event,a,nx,ne):
    global i
    j=i
    offset=10
    print("j:"+str(j)+"/"+str(ne))
    for j in range(0,int(ne),10): 
        plt.plot(a[j*nx:j*nx+nx])
        plt.draw()
    print(i)
    i+=offset
    if(i>ne):
        plt.close()
def main():
    
    file=h5py.File("./static_results.h5")
    d1=file.get("/numres")
    d2=file.get("/params1d")
    res=np.array(d1)
    param=np.array(d2)
    d1=0
    d2=0
    print(param)
    cev=functools.partial(drawsol,a=res,nx=int(param[4]),
ne=int(param[6]))
    plt.connect('button_press_event',cev)
    plt.show()

main()
