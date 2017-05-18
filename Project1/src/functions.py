import matplotlib.pyplot as plt
#%matplotlib inline
import math
import numpy as np
import sympy as sy
import pandas

def ff(x):
    return 100*math.exp(-10*x)
def exact(x):
    return 1.0-(1-math.exp(-10))*x-math.exp(-10*x)

i = 5
n = 10**i
h = 1.0/(n)
hh = h*h

d = np.zeros(n)
e = np.zeros(n)
f = np.zeros(n)
x = np.zeros(n)
u = np.zeros(n)

u[0] = u[n-1] = 0
for i in range(0,n,1):
    d[i] = 2

for i in range(0,n,1):
    x[i] = i*h
    f[i] = hh*ff(i*h)

for i in range(1,n,1):
    d[i] = (i+1)/i

for i in range(2,n,1):
    f[i] = f[i] + f[i-1]/d[i-1]

u[n-1] = f[n-1]/d[n-1]
for i in range(n-2,0,-1):
    u[i] = (f[i]+u[i+1])/d[i]

exact_results = []
for i in x:
    exact_results.append(exact(i))

def RelativeError(a,b):
    return math.log10(abs((a-b)/a))

errors = []
for i in range(1,n,1):
    errors.append(RelativeError(exact_results[i],u[i]))

lst1 = x[1:]
lst2 = u[1:]
lst3 = exact_results[1:]
lst4 = errors

table = pandas.DataFrame({'x': lst1, 'Computed': lst2, 'Exact': lst3,'Relative Error': lst4})

print(table)
