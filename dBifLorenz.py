import sys
import argparse
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def lorenz_attractor():

    steps = 10000
    dt = 0.01

    #sigma = 10
    rho = 28
    beta = 8.0/3 


    rho_min = 1#0.5
    rho_max = 13
    step = 1

    steps2 = (rho_max - rho_min)/step
    steps2 = int(steps2)
                        
    xs = np.empty((steps + 1,))
    ys = np.empty((steps + 1,))
    zs = np.empty((steps + 1,))
    
    #xs2 = np.empty((steps2 + 1))
    #ys2 = np.empty((steps2 + 1))
    #zs2 = np.empty((steps2 + 1))    
    xs2 = []
    ys2 = []
    zs2 = []
    
    print type(zs2)
    print np.shape(zs2)
    
    xs[0], ys[0], zs[0] = (0., 1., 1.05)

    # loop over steps in increments of dt
    
    k = 0
    #plt.ion()
    #plt.show()
    
    for sigma in np.arange(rho_min,rho_max,step):    
        print (rho)
        for i in np.arange(steps):
            xdot, ydot, zdot = lorenz(xs[i], ys[i], zs[i], rho, beta, sigma)
            xs[i+1] = xs[i] + (xdot * dt)
            ys[i+1] = ys[i] + (ydot * dt)
            zs[i+1] = zs[i] + (zdot * dt)              
        
        #print xs
        #print type(xs)
        #print np.shape(xs)
        #xs2[k] = xs
        #ys2[k] = ys
        #zs2[k] = zs
        #print(xs)
        
        xsc = xs.tolist()
        ysc = ys.tolist()
        zsc = zs.tolist()
        print(type(xsc))
        
        xs2.append(xsc)
        print('xs2')
        print(type(xs2))
        
        ys2.append(ysc)
        zs2.append(zsc)
        
        k = k+1
        
        #fig = plt.figure()
        #ax = Axes3D(fig)
        #ax.set_xlabel('X axis')
        #ax.set_ylabel('Y axis')
        #ax.set_zlabel('Z axis')
        #ax.plot(xs, ys, zs, label="lorenz")
        #ax.legend()
        #plt.draw()
        #print rho
        #print xs2
        #print type(xs2)
        #return xs2,ys2,zs2
        #if options.f:
	    #plt.savefig(options.f)
        #else:
        	#plt.show()
    #print np.shape(xs2)
    #print xs2[2]

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    str1 = (' sigma = from %.02f to%.02f , rho =  %0.02f , beta = %.02f , stepsize = %.02f' %(rho_min, rho_max,rho,beta,step))
    ax.set_title("Lorentz %s"%str1)
    colours=['r','g','m','b']
    #ax.legend()
    #plt.draw()
    colours2 = []
    
    range2use = len(xs2)/len(colours)
    print range2use
    for i in range(range2use):
        colours2.extend(colours)
    
    
    #ax.plot(xs2[5],ys2[5],zs2[5],colours2[5])
    
    for i in range(len(xs2)):
        ax.plot(xs2[i],ys2[i],zs2[i],colours2[i])
    
    #plt.draw()
    plt.show()
#plt.show()

def lorenz(x, y, z, r, b, s):
    xdot = s*(y - x)
    ydot = r*x - y - x*z
    zdot = x*y - b*z
    return xdot, ydot, zdot

def main():
    lorenz_attractor()

if __name__ == "__main__":
    main()
