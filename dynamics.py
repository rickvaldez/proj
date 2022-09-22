from pylab import *
import sys

t = arange(0,100) # time (generation number)
y = zeros(t.shape[0]) # vector for population size

mu = 3.5 # parameter value
y[0] = 0.75 # initial population 
#mu, y[0] = float(sys.argv[1]),float(sys.argv[2]) # parameters can be read 
# from the command line

for i in range(y.shape[0]-1): # evaluate logistic map for the desired number of generations
  y[i+1] = mu * y[i] * (1-y[i])


title('Growth parameter mu=%.2f and initial population X0=%.2f '%( mu,y[0]))
plot(t,y,'r') # plot population size vs time
xlabel('Population Generation n')
ylabel('Insect Population Xn')
show()
