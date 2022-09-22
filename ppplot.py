import matplotlib.pyplot as plt


workers = [1,2,3,4,5,6,7,8,9,10,11,12] 

time = [3.57223796844,1.79577302933,1.19228315353,1.03708791733,.894382953644,0.819687843323,0.979835033417,0.889533996582,0.875857830048,0.858896970749,0.815560102463,0.791584968567 ]


plt.plot(workers,time)
plt.title('Parallel Computing')
plt.xlabel('number of worker (# of cores)')
plt.ylabel('time seconds')
plt.show()
