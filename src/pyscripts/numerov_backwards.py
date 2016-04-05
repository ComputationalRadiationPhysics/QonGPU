import numpy as np
from matplotlib import pyplot as plt
import functools
it=0


def V(x):
    return 1/np.sqrt(1+x*x)


def Numerov(event):
    global it
    E=V(0)
    ne = 1e4
    dE=E/ne
    xm=10.0
    xmi=0.0
    nx = 1e4
    dx=(xm-xmi)/nx
    E=dE*it
    a=np.linspace(0,nx)
    a[-1]=0
    a[-2]=1e-10
    f1=0
    f2=0
    f3=0
    h=1/12*dx**2
    f=lambda x:E+V(x)
    for i in range(len(a)-1,-1,-1):
        f1=1+h*f(xm - dx*i)
        f2=1+h*f(xm - dx*(i-1))*5
        f3=1+h*f(xm - dx*(i-2))
        a[i-2]=1/f3*(2*f2*a[i-1]-f1*a[i])
        
    plt.plot(a)
    plt.draw()
    print("This is the: "+str(it)+"-th iteration")
    it+=1
def Numerov_forw(event):
    global it
    E=V(0)
    ne = 1e4
    dE=E/ne
    xm=10.0
    xmi=0.0
    nx = 1e4
    dx=(xm-xmi)/nx
    E=dE*it
    a=np.linspace(0,nx)
    a[0]=0
    a[1]=1e-10
    f1=0
    f2=0
    f3=0
    h=1/12*dx**2
    f=lambda x:E+V(x)
    for i in range(0,len(a)):    
        f1=1+h*f(xm - dx*i)
        f2=1+h*f(xm - dx*(i-1))*5
        f3=1+h*f(xm - dx*(i-2))
        a[i]=1/f1*(2*f2*a[i-1]-f3*a[i])
    print("Plotting...")
    plt.plot(a)
    plt.draw()
    print("This is the:"+str(it)+"iteration")
    it+=1
def main():
    xmax=10.0
    xmin=0.0
   
    plt.connect("button_press_event",Numerov)
    plt.show()
main()
