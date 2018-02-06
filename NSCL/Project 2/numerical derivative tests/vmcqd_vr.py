# -*- coding: utf-8 -*-
"""
Created on Thu May 25 12:45:51 2017

@author: ramir
"""
import numpy as np
import math
import time
from random import random
import sys
#from Solver import Solver
#
#from System import System

#from decimal import Decimal

def initialize(infile):
    '''Read line from file as parameters.
    The format:
        line 1: cycles, step length, number of particles, dimensions, w,
        starting value for alpha, starting value for beta, step size in alpha,
        step size in beta, number of alpha variations, number of beta variations
        line 2: jastrow factor enabled (bool), importance sampling enabled (bool)'''
    line_1 = infile.readline()
    params = line_1.split(',')
    cycles = eval(params[0])
    step_length = eval(params[1])
    n_particles = eval(params[2])
    dimensions = eval(params[3])
    w = eval(params[4])
    alpha_0 = eval(params[5])
    beta_0 = eval(params[6])
    alpha_step = eval(params[7])
    beta_step = eval(params[8])
    alpha_variations = eval(params[9])
    beta_variations = eval(params[10])
    line_2 = infile.readline()
    bool_params = line_2.split(',')
    jastrow_bool = eval(bool_params[0])
    importance_bool = eval(bool_params[1])
    #system = System(n_particles,dimensions,w, step_length)
    bool_params = [jastrow_bool, importance_bool]
    #hamiltonian_library = H(bool_params,system)
    #hamiltonian = hamiltonian_library.hamiltonians[bool_params]
    #solver = Solver(system, cycles,alpha_0,beta_0,alpha_step,beta_step,alpha_variations,beta_variations)
    return (cycles, step_length, n_particles, dimensions, w, alpha_0, beta_0, 
            alpha_step, beta_step, alpha_variations, beta_variations, 
            jastrow_bool, importance_bool)
    
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
    r_sum = 0.0    
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
    
#    harmonic oscillator  contribution
    for i in range(n_particles):
        r_single_particle = 0.0
        for j in range(dimensions):
            r_single_particle += r[i,j]**2
        e_potential += 0.5*r_single_particle

#    #Electron-electron contribution
#    for i1 in range(n_particles-1):
#        for i2 in range(i1+1,n_particles):
#            r_12 = 0.0
#            for j in range(dimensions):
#                r_12 += (r[i1,j] - r[i2,j])**2
#                
#            if jastrow_bool == True:
#                e_potential += 1/math.sqrt(r_12)
    
    return e_potential + e_kinetic

        
def get_params():
    '''Takes parameters from input file'''
    #infilename = input("Enter the input file name: ")
    args = sys.argv[1:]
    print(args)

    infilename= args[0]
    outfilename= args[1] 
    if infilename == '':
        print("No input file given in 1st argument")
    infile = open(infilename, 'r')
    print("Loading {}...".format(infilename))
    
    #outfilename = input("Enter the output file name: ")
    if outfilename == '':    
        print("No output file given in 2nd argument")
    outfile = open(outfilename, 'w')
    
    print("Writing on {}...".format(outfilename))
    print(infilename,outfilename)
    return infile, outfile
    

            

t_i = time.time()
#outfile = open(outfilename,'w')
h = .001
h2 = 1/(h**2)
alpha_step = 0.1
infile, outfile = get_params()

cycles, step_length, n_particles, dimensions, w, alpha_0, beta_0, alpha_step, beta_step, alpha_variations, beta_variations, jastrow_bool, importance_bool = initialize(infile)
#set up n x dim matrix representing positions
#r_ij represents the ith particle in the jth dimension
r_0 = np.zeros((n_particles,dimensions), np.double)
r_n = np.zeros((n_particles,dimensions), np.double)
#alphas = np.linspace(0.94,1.03,10)
#betas = np.linspace(0.36,0.45,10)
alpha_f = alpha_0+alpha_step*alpha_variations
beta_f = beta_0+beta_step*beta_variations
alphas = np.linspace(alpha_0,alpha_f,alpha_variations)
betas = np.linspace(beta_0,beta_step,beta_variations)
energy_min = 100.0

parameter_string = '''
***Parameters***
Monte Carlo cycles: {}
Step length: {}
{} Alpha variations in range:  {}-{}
{} Beta variations in range: {}-{}
Jastrow factor enabled? {}
Importance sampling enabled? {}
Number of particles: {}
Dimensions: {}
Oscillator Frequency w: {}
'''.format(cycles,step_length,alpha_variations,alpha_0,alpha_f,beta_variations,beta_0,beta_f,jastrow_bool,importance_bool,n_particles,dimensions,w)
print(parameter_string)

energy_matrix = np.zeros(shape=(alpha_variations,beta_variations))
print("Alpha, Beta, Energy, Variance, Error, Accept%")
for row in range(alpha_variations):
    alpha = alphas[row]
    for col in range(beta_variations):
        beta = betas[col]
        energy = energy2 = 0.0
        accept = 0.0
        delta_e = 0.0
        #print("Calculating alpha = {}, beta = {}...".format(alpha,beta))
        for i in range(n_particles):
            for j in range(dimensions):
                r_0[i,j] = step_length * (random() - .5)
                
        wf_0 = trial_wavefunction(w, alpha, r_0, n_particles, dimensions, beta, jastrow_bool)

        #Loop over cycles
        for cycle in range(cycles):        
            #Initial position

            #print(r_0)
    
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
        energy /= cycles
        energy2 /= cycles
        variance = abs(energy2 - energy**2)
            #print(energy, energy2, variance)
            #print(variance,type(variance))
        error = math.sqrt(variance/cycles);
            #...and write them to file
        #avg_energy = sum_energy/max_variations
        energy_matrix[row][col] = energy
        if energy < energy_min:
            energy_min = energy
            alpha_min = alpha
            beta_min = beta
        print('{:5} {:5} {:5} {:.5} {:.5} {:.5}\n'.format(alpha,beta,energy,variance,error,accept*1.0/(cycles)))
        outfile.write('%f %f %f %f %f %f\n' %(alpha,beta,energy,variance,error,accept*1.0/(cycles)))
        
outfile.close()
print("Minimum energy {} at alpha = {}, beta = {}".format(energy_min,alpha_min, beta_min))
print('\nDone. Results are in the outfile, formatted as:\n\
alpha, beta, <energy>, variance, error, acceptance ratio' )
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
