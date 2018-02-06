%matplotlib inline
import scipy.linalg
import numpy as np
import math
import time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    #function taken from
    #http://stackoverflow.com/questions/5842903/block-tridiagonal-matrix-python
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def plot(x, y):
    pass

#def time_take(n):
    #time_array = np.zeros(lim)
    #n_array = np.zeros(lim)
    a = [-1] * (n-1)
    b = [2] * n
    c = [-1] * (n-1)
    tic = time.clock()
    A = tridiag(a, b, c)
    B=scipy.linalg.lu(A)
    toc=time.clock()
    delta_t=toc-tic

    return delta_t

def LU_Decomp(n):
    a = [-1] * (n-1)
    b = [2] * n
    c = [-1] * (n-1)
    #tic = time.clock()
    A = tridiag(a, b, c)
    B=scipy.linalg.lu(A)
    #toc=time.clock()
    #delta_t=toc-tic

    return B

#change exp to 4 at your own caution.
#It takes absurdly long for reasons I cant figure out.
#The computation should be the same as above
exp = 3
trials = 20
n_array = np.zeros(exp-1)
#time_array = np.zeros(exp-1)

colors = ['r', 'b', 'g', 'y']
patches = []
for i in range(1,exp):
    #loop for within exponent
    n = 10**exp
    n_array[i-1] = i
    #time_array[i-1] = time_take(n)
    time_array = np.zeros(trials)
    for j in range(trials):
        #loop for each trial
        tic = time.clock()
        B = LU_Decomp(n)
        toc = time.clock()
        delta_t = toc-tic
        time_array[j] = delta_t

    print("10^" + str(i))
    patch = mpatches.Patch(color=colors[i-1], label="10^" + str(i))
    patches.append(patch)
    plt.scatter(np.arange(len(time_array)), time_array, color=colors[i-1])


plt.legend(handles=patches)

plt.title("Time (s) vs. N")
plt.ylabel("Time (s)")
plt.xlabel("Trials")
#print("The LU decomposition takes a time of :", sum(tim)/len(tim))
