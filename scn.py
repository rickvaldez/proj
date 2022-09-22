__author__ = 'rick'
from numpy import *
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as md
from scipy import stats


#x , t = loadtxt('AALK2214C038000_trades.csv',skiprows=1, usecols=(0,1), unpack = True)

size, px, timeStamp, bid, ask, stockPx, edge = genfromtxt('AALK2214C038000_trades.csv',skiprows=1,delimiter = ',', unpack = True)
#x.append(xt)

x = arange(938)
#print len(size)
#print size[937]
#print x
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

cent = ask - bid
#print cent

bid1 = (bid-min(bid))/ptp(bid)
ask1 = (ask-min(ask))/ptp(ask)
size1=(size-min(size))/ptp(size)#normalize
cent1 = (cent-min(cent))/ptp(cent)
px1 = (px - min(px))/ptp(px)
stockPx1 = (stockPx - min(stockPx))/ptp(stockPx)
#print size1
#normSize = size/(norm(size))

edge1 = edge * 5



ax.grid(b=True, color='0.75', which = 'both',axis = 'both',linestyle='-', linewidth =0.5)
#ax.plot(x,size,'-',color = 'm')
#ax.plot(x,bid,'o',color = 'b',s = 1 )
#ax.plot(x,ask,'o',color = 'r' , s = 1 )
ax.scatter(x,bid,color = 'b',s=2)
ax.scatter(x,ask,color = 'r',s=2)
#ax.plot(x,cent1,'-',color = 'm')
ax.plot(x,size1,'-',color = 'g')
ax.plot(x,edge1,'-',color = 'c')
ax.plot(x,px,'-',color = 'y')
ax.plot(x,stockPx1,'-',color = 'k')
plt.show()
