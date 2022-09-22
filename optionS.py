import lxml.html
import calendar 
from datetime import datetime
from numpy import array

exDate = "2015-04-10"
symbol = "aapl"
dt = datetime.strptime(exDate,'%Y-%m-%d')
ym = calendar.timegm(dt.utctimetuple())
#print ym
url = 'http://finance.yahoo.com/q/op?s=%s&date=%s' % (symbol,ym,)
doc = lxml.html.parse(url)
table =doc.xpath('//table[@class="details-table quote-table Fz-m"]/tbody/tr')

rows = []
uoa = []

for tr in table:
    d = [td.text_content().strip().replace(',','') for td in tr.xpath('./td')]
    rows.append(d)


print rows[22][8]



#print len(rows)

l = len(rows)

"""
for i in range (0,l):

	openInt = rows[i][8]
	vol = rows[i][7]
	
	print '\n ....vol', vol
	print 'openint', openInt
		
	if float(vol) > float(openInt)*2.5:		
		print vol, openInt
		uoa.append(i)
print uoa
"""


for i in range (0,l):

	openInt = rows[i][8]
	vol = rows[i][7]
	
	print '\n ....vol', vol
	print 'openint', openInt
		
	if float(vol) > float(openInt)*2.5:		
		print vol, openInt
		uoa.append(i)
print uoa



for j in uoa:
	print rows[j]#[7:9]





#openInt = rows[22][8]
#vol = rows[22][7]

#print openInt
#print vol



