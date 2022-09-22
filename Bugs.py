"""
Updated bug population bifurcation diagram script
"""  

from pylab import *

#ion()

# parameter m is in [m_min,m_max]
m_min = 1.0
m_max = 4.0
step = 0.0025 # step size in m for the diagram

y = 0.5 # initial population
lasty = int(1000 * y) # to avoid plotting a point twice (at least in a row)
N = 5 # a point will be put on the diagram every N iterations
count=0
for m in arange(m_min, m_max, step):
    for i in range(400): # to avoid transients
        y = m*y*(1-y) # applying logistic map   
    for i in range(200): # next 200 generations would go to the diagram
        y = m*y*(1-y)   
        inty = int(1000 * y) # keep three significant digits
        if inty != lasty and count%N==0: # no repeats, check generation number
            plot(m,y,'ro',markersize = 1,markerfacecolor = 'r',markeredgecolor = 'r') #'g.') # put a dot on the diagram
        lasty = inty
        count+=1
title('Bifurcation plot of attractor population X* vs growth rate mu')
xlabel('mu')
ylabel('x*')
show()
