# -*- coding: utf-8 -*-
"""
Created on Thu May 25 12:45:51 2017

@author: ramir
"""
import numpy as np
import math
import time
from random import random
#from Solver import Solver
#
#from System import System

#from decimal import Decimal

def initialize(infilename):
    '''Read line from file as parameters.
    The format is n_particles, dimensions, max_variations, MC cycles, step length'''
    with open(infilename,'r') as file:
        line = file.readline()
        params = line.split(',')
        max_variations = eval(params[0])
        n_cycles = eval(params[1])
        step_length = eval(params[2])
        alpha = eval(params[3])
        beta = eval(params[4])
        jastrow_bool = eval(params[5])
        n_particles = eval(params[6])
        dimensions = eval(params[7])
        w = eval(params[8])
        #system = System(n_particles,dimensions,w)
        #solver = Solver(variations,step_length)
        
    #return system, solver
    return max_variations, n_cycles, step_length,alpha,beta,jastrow_bool,n_particles,dimensions,w
    
def jastrow(wf,a,r,beta,n_particles,dimension):
    #jastrow_factor = 1.0
    for i in range(n_particles-1):
        for j in range(i+1,n_particles):
            r_12 = 0.0
            for k in range(dimension):
                r_12 += (r[i,k] - r[j,k])**2
            arg = math.sqrt(r_12)
            wf *= math.exp(arg/(1.0+beta*arg))
    return wf
    
    
def trial_wavefunction(w, alpha, r, n_particles, dimensions, beta, jastrow_bool=False):
    r_sum = 0    
    for i in range(n_particles):
        #loop over each particle
        r_ij_particle = 0.0
        for j in range(dimensions):
            r_ij_particle += r[i,j]**2
        r_sum += r_ij_particle
    wf = math.exp(-0.5*alpha*w*r_sum)
    if jastrow_bool == True:
        wf = jastrow(wf,alpha, r, beta, n_particles, dimensions)
    return wf
        
        
def local_energy(w, alpha, r,n_particles,dimensions, wf, beta, jastrow_bool, h=0.001, h2 = 1E6):
    #Kinetic energy
    r_plus = r.copy()
    r_minus = r.copy()
    e_kinetic = 0.0
    for i in range(n_particles):
        for j in range(dimensions):
            r_plus[i,j] = r[i,j] + h
            r_minus[i,j] = r[i,j] - h
            wf_minus = trial_wavefunction(w, alpha, r_minus, dimensions, n_particles, beta, jastrow_bool)
            wf_plus = trial_wavefunction(w, alpha, r_plus, dimensions, n_particles, beta, jastrow_bool)
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
                
            if jastrow_bool == True:
                e_potential += 1/math.sqrt(r_12)
    
    return e_potential + e_kinetic
    
#def optimize_parameters(alphas, betas):
#    for alpha in alphas:
#        for beta in betas:
#            sum_energy = 0.0
#            print("Calculating alpha = {}, beta = {}".format(alpha,beta))
        
def get_params():
    infilename = input("Enter the input file name: ")
    if infilename == '':
        infilename = 'in.txt'
    
    outfilename = input("Enter the output file name: ")
    if outfilename == '':    
        outfilename = 'out.txt'
        
    jastrow_str = input("Jastrow? [y/n]: ")
    if jastrow_str == 'y':
        jastrow_bool = True
    else:
        jastrow_bool = False
        
    return infilename, outfilename, jastrow_bool
    
infilename, outfilename, jastrow_bool = get_params()
        

t_i = time.time()
outfile = open(outfilename,'w')
h = .001
h2 = 1/(h**2)
alpha_step = 0.1
max_variations, n_cycles, step_length,alpha,beta,jastrow_bool,n_particles,dimensions,w = initialize(infilename)
#set up n x dim matrix representing positions
#r_ij represents the ith particle in the jth dimension
r_0 = np.zeros((n_particles,dimensions), np.double)
r_n = np.zeros((n_particles,dimensions), np.double)
#alphas = np.linspace(0.94,1.03,10)
#betas = np.linspace(0.36,0.45,10)
alphas = np.linspace(0.95,1.05,11)
betas = np.linspace(0.35,0.45,11)
energy_min = 100.0
print("Parameters: \nnumber of particles = {} \ndimensions = {} \nvariations = {} \nMC cycles = {} \nstep length = {} \nw = {} \nJastrow Factor enabled? {}".format(n_particles, dimensions, max_variations, n_cycles, step_length, w, jastrow_bool))
for alpha in alphas:
    for beta in betas:
        sum_energy = 0.0
        print("Calculating alpha = {}, beta = {}...".format(alpha,beta))
        for variation in range(max_variations):
            #alpha += alpha_step
        
            energy = energy2 = 0.0
            accept = 0.0
            delta_e = 0.0
        
            #Initial position
            for i in range(n_particles):
                for j in range(dimensions):
                    r_0[i,j] = step_length * (random() - .5)
            #print(r_0)
            wf_0 = trial_wavefunction(w, alpha, r_0, n_particles, dimensions, beta, jastrow_bool)
            
            #Loop over MC cycles
            for cycle in range(n_cycles):
        
                #Trial position
                for i in range(n_particles):
                    for j in range(dimensions):
                        r_n[i,j] = r_0[i,j] + step_length * (random() - .5)
        
                wf_n = trial_wavefunction(w, alpha, r_n, n_particles, dimensions, beta, jastrow_bool)
                
                #Metropolis test to see whether we accept the move
                if random() < wf_n**2 / wf_0**2:
                    r_0 = r_n.copy()
                    wf_0 = wf_n
                    accept += 1
        #            
                #update expectation values
                delta_e = local_energy(w,alpha,r_0,n_particles,dimensions,wf_0, beta, jastrow_bool)
                energy += delta_e
                energy2 += delta_e**2
        #
            #We calculate mean, variance and error ...
            energy /= n_cycles
            sum_energy += energy
            energy2 /= n_cycles
            variance = abs(energy2 - energy**2)
            #print(energy, energy2, variance)
            #print(variance,type(variance))
            error = math.sqrt(variance/n_cycles);
            #...and write them to file
        avg_energy = sum_energy/max_variations
        if avg_energy < energy_min:
            energy_min = avg_energy
            alpha_min = alpha
            beta_min = beta
        outfile.write('%f %f %f %f %f %f\n' %(alpha,beta,avg_energy,variance,error,accept*1.0/(n_cycles)))
        
outfile.close()
print("Minimum energy {} at alpha = {}, beta = {}".format(energy_min,alpha_min, beta_min))
print('\nDone. Results are in the file "%s", formatted as:\n\
alpha, beta, <energy>, variance, error, acceptance ratio' %(outfilename))
t_f = time.time()
delta_t = t_f - t_i
print("computation time took {}s".format(delta_t))

#jastrow
#alpha ~ 0.98, beta ~ 0.4, energy close to 3
#take away jastrow, calculate energy
#importance sampling
#calculation of covariance, standard deviation with blocking
#finding minimum fr multiple functions, find optimal alpha and beta
#conjugate gradient
#parallelize
