# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:33:50 2017

@author: ramir
"""

import numpy as np
import matplotlib.pyplot as plt

file = open("out.txt","r")
alpha = []
pots = []
variances = []
errors = []
acceptance= []
for line in file:
    values = line.split()
    alpha.append(float(values[0]))
    pots.append(float(values[1]))
    variances.append(float(values[2]))
    errors.append(float(values[3]))
    acceptance.append(float(values[4]))

    
#plt.plot(alpha,pots)
plt.errorbar(alpha, pots, yerr=errors,ecolor='r')
plt.xlabel("alpha")
plt.ylabel("potential")
plt.show()