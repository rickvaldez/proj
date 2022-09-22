
 #stability test for nu > .25
# Previous value nu = 0.09 
 
#Old kappa 237 , new kappa 650 
# 
 
from numpy import *
import matplotlib.pylab as p
from mpl_toolkits.mplot3d import Axes3D 

Nx = 101;        Nt = 15000;     Dx = 0.03;     Dt = 0.9                                                              
KAPPA = 1250.; SPH = 900.; RHO = 2700. # conductivity, specf heat, density                                                      
T = zeros( (Nx, 2), float);  Tpl = zeros( (Nx, 31), float)  
                                     
print("Working, wait for figure after count to 10")

for ix in range (1, Nx - 1):  T[ix, 0] = 100.0;     # initial temperature
T[0,0] = 0.0 ;   T[0,1] = 0.                 # first and last points at 0
T[Nx-1,0] = 0. ; T[Nx-1,1] = 0.0
cons = KAPPA/(SPH*RHO)*Dt/(Dx*Dx);                             # constant
m = 1                        # counter for rows, one every 300 time steps

for t in range (1, Nt):                                  # time iteration
   for ix in range (1, Nx - 1):                      # Finite differences
      T[ix, 1] = T[ix, 0] +  cons*(T[ix+1, 0] + T[ix-1, 0] - 2.*T[ix,0])                                                        
   if t%500 == 0 or t == 1:          # for t = 1 and every 300 time steps
        for ix in range (1, Nx - 1, 2): Tpl[ix, m] = T[ix, 1]   
        print(m)   
        m = m + 1                       # increase m every 300 time steps
   for ix in range (1, Nx - 1):  T[ix, 0] = T[ix, 1]# 100 positons at t=m
x = list(range(1, Nx - 1, 2))                  # plot every other x point
y = list(range(1, 30))                      # every 10 points in y (time)
X, Y = p.meshgrid(x, y)                      # grid for position and time

def functz(Tpl):                           # Function returns temperature
    z = Tpl[X, Y]       
    return z




Z = functz(Tpl)              
fig = p.figure()                                          # create figure
ax = Axes3D(fig)                                             # plots axis
ax.plot_wireframe(X, Y, Z, color = 'r')                   # red wireframe
ax.set_xlabel('Position')                                    # label axes
ax.set_ylabel('time')
ax.set_zlabel('Temperature')
p.show()                               # shows figure, close Python shell
print("finished")                               
