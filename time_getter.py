
import datetime
import sys
import numpy as np

fname = sys.argv[1]

if fname[0]=='.'
    print('File: ', fname)
    print('Algorithm will fail. Do not use a leading dot on the filename.')

o = fname.split('.')[0][-6:]
h = o[:2]
m = o[2:4]
s = o[4:]
print(h,m,s)

s = float(s) + (float(h)*60+float(m))*60

# s is now the base time to be added to all timestamps in the data.

