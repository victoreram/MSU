# -*- coding: utf-8 -*-
"""
Created on Thu May 25 12:45:51 2017

@author: ramir
"""
import numpy as np
import math
import time
from random import random
#from decimal import Decimal

def init(input_file):
    '''Read line from file as parameters.
    The format is n_particles, dimensions, max_variations, MC cycles, step length'''
    with open(input_file,'r') as file:
        line = file.readline()
        params = line.split(',')
        n_particles = eval(params[0])
        dimensions = eval(params[1])
        max_variations = eval(params[2])
        n_cycles = eval(params[3])
        step_length = eval(params[4])
        alpha = eval(params[5])
        alpha_step = eval(params[6])
        w = eval(params[7])
    return n_particles, dimensions, max_variations, n_cycles, step_length, alpha, alpha_step,w
    
def jastrow(wf,a,r,beta,n,dimension):
    for i in range(n-1):
        for j in range(i+1,n):
            r_12 = 0.0
            for k in range(dimension):
                r_12 += (r[i,k] - r[j,k])**2
            r_12 = math.sqrt(r_12)
            wf *= math.exp(.5*r_12)
    return wf
    
def energy_l1(w,r_1,r_2,alpha):
    return .5*w**2*(r_1**2+r_2**2)*(1-alpha**2) + 2*alpha*w
    
def trial_wavefunction(w, alpha, r, n_particles, dimensions, beta = 0.3, jastrow_factor=False):
    r_sum = 0    
    for i in range(n_particles):
        #loop over each particle
        r_ij_particle = 0.0
        for j in range(dimensions):
            r_ij_particle += r[i,j]**2
        r_sum += r_ij_particle
    wf = math.exp(-0.5*alpha*w*r_sum)
    if jastrow_factor == True:
        wf *= jastrow(wf,alpha, r, beta, n_particles, dimensions)
    return wf
        
        
def local_energy(w, alpha, r,n_particles,dimensions, wf, h=0.001, h2 = 1E6):
    #Kinetic energy
    r_plus = r.copy()
    r_minus = r.copy()
    e_kinetic = 0.0
    for i in range(n_particles):
        for j in range(dimensions):
            r_plus[i,j] = r[i,j] + h
            r_minus[i,j] = r[i,j] - h
            wf_minus = trial_wavefunction(w, alpha, r_minus, dimensions, n_particles)
            wf_plus = trial_wavefunction(w, alpha, r_plus, dimensions, n_particles)
            e_kinetic -= wf_minus+wf_plus-2*wf;
            r_plus[i,j] = r[i,j]
            r_minus[i,j] = r[i,j]
    
    e_kinetic = .5*h2*e_kinetic/wf
    #print("e_kinetic",  e_kinetic)
    #Potential energy
    e_potential = 0.0
    
    #harmonic oscillator  contribution
    for i in range(n_particles):
        r_single_particle = 0.0
        for j in range(dimensions):
            r_single_particle += r[i,j]**2
        e_potential += 0.5*r_single_particle

    #Electron-electron contribution
    for i1 in range(n_particles-1):
        for i2 in range(i1+1,n_particles):
            r_12 = 0.0
            for j in range(dimensions):
                r_12 += (r[i1,j] - r[i2,j])**2
                
#            e_potential += 1/math.sqrt(r_12)
    
    return e_potential + e_kinetic
input_file = input("Enter the input file name: ")
if input_file == '':
    input_file = 'infile.txt'

outfilename = input("Enter the output file name: ")
if outfilename == '':    
    outfilename = 'outfile.txt'
t_i = time.time()
outfile = open(outfilename,'w')
h = .001
h2 = 1/(h**2)
n_particles, dimensions, max_variations, n_cycles, step_length, alpha, alpha_step,w = init(input_file)
#set up n x dim matrix representing positions
#r_ij represents the ith particle in the jth dimension
r_0 = np.zeros((n_particles,dimensions), np.double)
r_n = np.zeros((n_particles,dimensions), np.double)
for variation in range(max_variations):
    alpha += alpha_step

    energy = energy2 = 0.0
    accept = 0.0
    delta_e = 0.0

    #Initial position
    for i in range(n_particles):
        for j in range(dimensions):
            r_0[i,j] = step_length * (random() - .5)
    #print(r_0)
    wf_0 = trial_wavefunction(w, alpha, r_0, n_particles, dimensions)
    
    
    #Loop over MC cycles
    for cycle in range(n_cycles):

        #Trial position
        for i in range(n_particles):
            for j in range(dimensions):
                r_n[i,j] = r_0[i,j] + step_length * (random() - .5)

        wf_n = trial_wavefunction(w, alpha, r_n, n_particles, dimensions)
        
        #Metropolis test to see whether we accept the move
        if random() < wf_n**2 / wf_0**2:
            r_0 = r_n.copy()
            wf_0 = wf_n
            accept += 1
#            
        #update expectation values
        delta_e = local_energy(w,alpha,r_0,n_particles,dimensions,wf_0)
        energy += delta_e
        energy2 += delta_e**2
#
    #We calculate mean, variance and error ...
    energy /= n_cycles
    energy2 /= n_cycles
    variance = energy2 - energy**2
    #variance_scalar = np.mean(variance)
    #energy_scalar = np.mean(energy)
    #print("<Energy> {}, <Energy>^2: {}, Variance: {} ".format(energy, energy2, variance_scalar))
    #print(variance, sum(variance))
    error = np.sqrt(variance/n_cycles)
    #print("error", error)
    #...and write them to file
    outfile.write('%f %f %f %f %f\n' %(alpha,energy,variance,error,accept*1.0/(n_cycles)))

outfile.close()

print('\nDone. Results are in the file "%s", formatted as:\n\
alpha, <energy>, variance, error, acceptance ratio' %(outfilename))
t_f = time.time()
delta_t = t_f - t_i
print("computation time took {}s".format(delta_t))

#jastrow
#importance sampling
#calculation of covariance, standard deviation with blocking
#finding minimum fr multiple functions, find optimal alpha and beta
#conjugate gradient
#parallelize
