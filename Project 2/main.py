# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 14:31:37 2017

@author: ramir
"""
#import numpy as np
#import math
#from random import random
from Solver import Solver
from System import System
#from Hamiltonians import Hamiltonian as H
import time


def initialize(infile):
    '''Read line from file as parameters.
    The format is cycles, step length, alpha, beta, 
    jastrow factor (True or False), number of particles, dimensions, 
    oscillator frequency w'''
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
    system = System(n_particles,dimensions,w, step_length)
    solver = Solver(cycles,alpha_0,beta_0,alpha_step,beta_step,alpha_variations,beta_variations,jastrow_bool)
    return solver, system

def get_params():
    infilename = input("Enter the input file name: ")
    if infilename == '':
        infilename = 'in.txt'
    infile = open(infilename, 'r')
    print("Loading {}...".format(infilename))
    
    outfilename = input("Enter the output file name: ")
    if outfilename == '':    
        outfilename = 'out.txt'
    outfile = open(outfilename, 'w')
    
    print("Writing on {}...".format(outfilename))
    return infile, outfile
    
infile, outfile = get_params()
solver, system = initialize(infile)
parameter_string = '''
***Parameters***
Monte Carlo cycles: {}
Step length: {}
{} Alpha variations in range:  {}-{}
{} Beta variations in range: {}-{}
Jastrow factor enabled? {}
Number of particles: {}
Dimensions: {}
Oscillator Frequency w: {}
'''.format(solver.mc_cycles,system.step_length,solver.alpha_variations,solver.alpha,solver.final_alpha,solver.beta_variations,solver.beta,solver.final_beta,solver.jastrow_bool,system.number_of_particles,system.dimensions,system.w)
print(parameter_string)
t_i = time.time()

solver.optimize_parameters(system, outfile)
infile.close()
outfile.close()
print('\nDone. Results are in the output file, formatted as:\n\
alpha, beta, <energy>, variance, error, acceptance ratio' )

t_f = time.time()
delta_t = t_f - t_i
print("computation time took {}s".format(delta_t))
#take away jastrow, calculate energy
#importance sampling
#calculation of covariance, standard deviation with blocking
#finding minimum fr multiple functions, find optimal alpha and beta
#conjugate gradient
#parallelize

