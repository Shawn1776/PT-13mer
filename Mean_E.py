import re # for rstrip()

f=open("./T_1_STEP_10E7.txt","r")
tt = f.readlines()
f.close()
#len(tt)
t_len = len(tt)/2   # drop the first half of data
t_len = int(t_len)  # cast to int type

t_len = int(100000)
t_rest = tt[t_len:]

#len(t_rest)

t_removeNL=[]
t_removeNL=[x.rstrip() for x in t_rest] # remove the "\n" for each line
#len(t_removeNL)

 
tu_int = list(map(lambda x:(float(x.split("\t")[0]),float(x.split("\t")[1])), t_removeNL)) #
step = [x[0] for x in tu_int]
E = [x[1] for x in tu_int] # only E, use x[] not x()

E_avg = sum(E) / float(len(E))

print("E_avg is: %s" % E_avg)

import pylab as pl

pl.plot(step,E)

pl.show()
