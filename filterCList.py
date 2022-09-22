__author__ = 'rick'


import csv


look_for = set(['41' , '42' , '57' , '37' , '38' , '56'])

with open('expired-subs-2014-09-28.csv','rb') as inf, open('filt-expired-subs-2014-09-28.csv','wb') as outf:
    incsv = csv.reader(inf, delimiter=',',skipinitialspace=True)
    outcsv = csv.writer(outf, delimiter=',')
    outcsv.writerows(row for row in incsv if row[6] in look_for)
