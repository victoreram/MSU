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


def initialize(infile):
    '''Read line from file as parameters.
    The format is variations, cycles, step length, alpha, beta, 
    jastrow factor (True or False), number of particles, dimensions, 
    oscillator frequency w'''
    line = infile.readline()
    params = line.split(',')
    variations = eval(params[0])
    cycles = eval(params[1])
    step_length = eval(params[2])
    alpha_0 = eval(params[3])
    beta_0 = eval(params[4])
    jastrow_bool = eval(params[5])
    n_particles = eval(params[6])
    dimensions = eval(params[7])
    w = eval(params[8])
    alpha_step = eval(params[9])
    beta_step = eval(params[10])
    beta_variations = eval(params[11])
    system = System(n_particles,dimensions,w)
    solver = Solver(variations,cycles,step_length,alpha_0,beta_0,alpha_step,beta_step,beta_variations,jastrow_bool)
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
    
    print("Writing on {}".format(outfilename))
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
'''.format(solver.mc_cycles,solver.step_length,solver.alpha_variations,solver.alpha,solver.final_alpha,solver.beta_variations,solver.beta,solver.final_beta,solver.jastrow_bool,system.number_of_particles,system.dimensions,system.w)
print(parameter_string)
solver.optimize_parameters(system, outfile)
infile.close()
outfile.close()
print('\nDone. Results are in the output file, formatted as:\n\
alpha, beta, <energy>, variance, error, acceptance ratio' )
#take away jastrow, calculate energy
#importance sampling
#calculation of covariance, standard deviation with blocking
#finding minimum fr multiple functions, find optimal alpha and beta
#conjugate gradient
#parallelize

