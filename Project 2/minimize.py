# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:16:55 2017

@author: ramir
"""
import numpy as np
infilename = "outjas.txt"
infile = open(infilename, "r")
alphas = np.linspace(0.9,1.1,21)
betas = np.linspace(0.3,0.5,21)
minima = 5.0
for line in infile:
    #print(line)
    params = line.split()
    alpha = float(params[0])
    beta = float(params[1])
    energy = float(params[2])
    if energy < minima:
        minima = energy
        alpha_min = alpha
        beta_min = beta

print(minima, alpha_min, beta_min)
#print(energies)