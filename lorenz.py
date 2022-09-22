import sys
import argparse
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def lorenz_attractor():

    steps = 10000
    dt = 0.01

    sigma = 10
    rho = 16
    beta = 8.0/3 


                        
    xs = np.empty((steps + 1,))
    ys = np.empty((steps + 1,))
    zs = np.empty((steps + 1,))
    xs[0], ys[0], zs[0] = (0,1,1.05)#(0., 1., 1.05)

    # loop over steps in increments of dt
    for i in np.arange(steps):
        xdot, ydot, zdot = lorenz(xs[i], ys[i], zs[i], rho, beta, sigma)
        xs[i+1] = xs[i] + (xdot * dt)
        ys[i+1] = ys[i] + (ydot * dt)
        zs[i+1] = zs[i] + (zdot * dt)              

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.plot(xs, ys, zs, label="lorenz")

    str1 = (' sigma = %.02f , rho = %.02f , beta = %.02f' %(sigma,rho,beta))
    ax.set_title("Lorentz %s"%str1)
#    ax.text2D(0.05,0.095,"info is %s"%str1,transform = ax.transAxes)
    ax.legend()
    plt.show()

def lorenz(x, y, z, r, b, s):
    xdot = s*(y - x)
    ydot = r*x - y - x*z
    zdot = x*y - b*z
    return xdot, ydot, zdot

def main():
    lorenz_attractor()

if __name__ == "__main__":

    main()
