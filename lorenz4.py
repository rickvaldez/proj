
import numpy as numpy
import matplotlib.pyplot 

def lorenz():
    h = 0.001
    end_time = 50.
    num_steps = int(end_time / h)

    sigma = 10.
    beta = 8. / 3.
    rho = 28

    x = numpy.zeros([num_steps + 1, 2])
    y = numpy.zeros([num_steps + 1, 2])
    z = numpy.zeros([num_steps + 1, 2])
    times = h * numpy.array(range(num_steps + 1))    

    x[0, 0] = 0.
    y[0, 0] = 0.3 
    z[0, 0] = 40.

    x[0, 1] = 0.
    y[0, 1] = 0.300000000000001 
    z[0, 1] = 40.

    distance = numpy.zeros(num_steps + 1)
    for step in range(num_steps):
        x[step+1, 0] = x[step, 0] + h*sigma*(y[step, 0] - x[step, 0])
        y[step+1, 0] = y[step, 0] + h*(x[step, 0]*(rho-z[step, 0])-y[step, 0])
        z[step+1, 0] = z[step, 0] + h*(x[step, 0]*y[step, 0] - beta*z[step, 0])
        x[step+1, 1] = x[step, 1] + h*sigma*(y[step, 1] - x[step, 1])
        y[step+1, 1] = y[step, 1] + h*(x[step, 1]*(rho-z[step, 1])-y[step, 1])
        z[step+1, 1] = z[step, 1] + h*(x[step, 1]*y[step, 1] - beta*z[step, 1])
        
        distance[step+1] = abs(x[step+1,0] - x[step+1, 1])
        
    return times, distance, x, z, y , sigma , rho

times, distance, x, z, y, sigma ,rho = lorenz()

#show_plot
def plot_me():
    axes = matplotlib.pyplot.gca()

    matplotlib.pyplot.plot(times, x)
    str1 = ('sigma = %.03f rho = %.03f , beta = 8/3' %(sigma,rho)) 
    matplotlib.pyplot.title('Transient for %s'%(str1))
    axes.set_xlabel('t')
    axes.set_ylabel('Distance of x values')
    
    # 
    #matplotlib.pyplot.plot(x[:, 0], z[:, 0])
    #matplotlib.pyplot.plot(x[:, 1], z[:, 1])
    #axes.set_xlabel('x')
    #axes.set_ylabel('z')
    matplotlib.pyplot.show()

plot_me()
