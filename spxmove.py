__author__ = 'rick'

import numpy as np

data = np.genfromtxt('spx.csv', delimiter = ',', skiprows = 1, usecols = 4)

dateData = np.genfromtxt('spx.csv', delimiter = ',', skiprows = 1, usecols = 0, dtype = str)
#print data

percentDiff = np.diff(data) / data[:-1] * 100

np.set_printoptions(threshold='nan')

percentChange = []
indexForDate = []

print percentDiff, ' percent diff'
j= 0;
for x in percentDiff:
    if abs(x) > 1 :
        percentChange.append(x)

        indexForDate.append(j)
        j = j + 1
    else:
        j = j + 1


#print percentChange

print indexForDate

consecutiveDates = zip(indexForDate,indexForDate[1:],indexForDate[2:])



#[(a,b) for a,b in zip(indexForDate,indexForDate[1:]) if b == a+1:]

print consecutiveDates


consecutiveDatesInRow = []

for a, b, c in consecutiveDates:
        if c == b +1 and b == a +1:
            appendVal = a
            consecutiveDatesInRow.append(appendVal)

print consecutiveDatesInRow

#print d

for y in consecutiveDatesInRow:
    print dateData[y]
