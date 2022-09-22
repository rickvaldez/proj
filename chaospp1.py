import numpy
import pylab
import sys
import pp





def chaos(m,y,lasty,N,count):

    
    popu = []
    for i in range(400):
	    y = m*y*(1-y)
    for i in range(200):
    	y = m*y*(1-y)
    	inty = int(1000*y)
        if inty != lasty and count%N==0:
            popu.append([m,y])
        lasty = inty
        count+=1
    return (popu)

mValuesPassed = []
for i in numpy.arange(1,4,0.0025):
    mValuesPassed.append(i)


y = 0.5
lasty = int(1000*y)
N = 5
count = 0

######### PP part

print """Usage: python sum_primes.py [ncpus]
    [ncpus] - the number of workers to run in parallel,
    if omitted it will be set to the number of processors in the system"""

# tuple of all parallel python servers to connect with
#ppservers = ("curie:35000","bohr:35002","pauli:35002","rabi:35002")
#ppservers = ("127.0.0.1:60000", )
ppservers = ()

if len(sys.argv) > 1:
    ncpus = int(sys.argv[1])
    # Creates jobserver with ncpus workers
    job_server = pp.Server(ncpus, ppservers=ppservers, secret="blah")
else:
    # Creates jobserver with automatically detected number of workers
    job_server = pp.Server(ppservers=ppservers, secret="blah")

print "Starting pp with", job_server.get_ncpus(), "workers"


########

points = []




job1 = [(inp,job_server.submit(chaos,(inp,y,lasty,N,count, ),(),())) for inp in mValuesPassed]



print job1[25] 

for inp, job in job1:
    points.extend(job())
    #print 'inp',inp,'job',job()
  
print (len(points))
#print result

'''
for i in mValuesPassed: 
    popu = chaos(i,y,lasty,N,count)
    points.extend(popu)
    
print(len(points))
print(points[50][0],points[50][1])
'''

integ = 0

for i in points:
    pylab.plot(points[integ][0],points[integ][1],'g.')
    integ = integ+1

pylab.show()

job_server.print_stats()


