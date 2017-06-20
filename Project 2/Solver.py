# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:42:20 2017

@author: ramir
"""
import math
import numpy as np
#from System import System
from random import random
class Solver:
    def __init__(self, alpha_variations, mc_cycles, step_length, alpha, beta, alpha_step,beta_step,beta_variations,jastrow_bool):
        '''Solver(variations, mc_cycles, step_length, alpha, beta, jastrow factor enabled
        (True or False))'''
        self.alpha_variations = alpha_variations
        self.beta_variations = beta_variations
        self.mc_cycles = mc_cycles
        self.step_length = step_length
        self.alpha = alpha
        self.beta = beta
        self.jastrow_bool = jastrow_bool
        self.alpha_step = alpha_step
        self.beta_step = beta_step
        self.final_alpha = alpha + (alpha_variations-1)*alpha_step
        self.final_beta = beta + (beta_variations-1)*beta_step
        self.alphas = np.linspace(alpha,self.final_alpha,alpha_variations)
        self.betas = np.linspace(beta,self.final_beta,beta_variations)
        
        
    def optimize_parameters(self, system, outfile, minimum_energy=1000.0):
        alphas = self.alphas.copy()
        betas = self.betas.copy()
        for alpha in alphas:
            for beta in betas:
                #find minimum
                self.alpha = alpha
                self.beta = beta
                print("calculating alpha = {}, beta = {}...".format(alpha,beta))
                results = self.MC_calculations(system)
                self.write_to_file(outfile,results)
                energy = results[2]
                if energy < minimum_energy:
                    minimum_energy = energy
                    minimum_alpha = results[0]
                    minimum_beta = results[1]
        print("Minimum energy {} at alpha = {}, beta = {}".format(minimum_energy, minimum_alpha, minimum_beta))
        return minimum_alpha, minimum_beta, minimum_energy                      
            
    def trial_wavefunction(self, r, system):
        r_sum = 0.0
        for i in range(system.number_of_particles):
            #loop over each particle
            r_ij_particle = 0.0
            for j in range(system.dimensions):
                r_ij_particle += r[i,j]**2
            r_sum += r_ij_particle
        wf = math.exp(-0.5*self.alpha*system.w*r_sum)
        if self.jastrow_bool == True:
            wf = self.jastrow(wf, r, system)
        return wf
        
    def jastrow(self,wf,r,system):
        for i in range(system.number_of_particles-1):
            for j in range(i+1,system.number_of_particles):
                r_12 = 0.0
                for k in range(system.dimensions):
                    r_12 += (r[i,k] - r[j,k])**2
                arg = math.sqrt(r_12)
                wf *= math.exp(0.5*arg/(1.0+self.beta*arg))
        return wf
        #analytic expression version
        #for i in range(system.number_of_particles)
        
    def local_energy(self,r, wf, system, h=0.001, h2 = 1E6):
        #Kinetic energy
        n_particles = system.number_of_particles
        dimensions = system.dimensions
        r_plus = r.copy()
        r_minus = r.copy()
        e_kinetic = 0.0
        for i in range(n_particles):
            for j in range(dimensions):
                r_plus[i,j] = r[i,j] + h
                r_minus[i,j] = r[i,j] - h
                wf_minus = self.trial_wavefunction(r_minus, system)
                wf_plus = self.trial_wavefunction(r_plus, system)
                e_kinetic -= wf_minus+wf_plus-2*wf;
                r_plus[i,j] = r[i,j]
                r_minus[i,j] = r[i,j]
        
        e_kinetic = .5*h2*e_kinetic/wf
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
                if self.jastrow_bool == True:
                    e_potential += 1.0/math.sqrt(r_12)
        
        return e_potential + e_kinetic
        
    def local_energy1(self, system):
        w = system.w
        energy = 0.0
        r_squared = 0.0
        for particle in system.particles:
            r_i = particle.r_squared()
            r_squared += r_i
            
        energy = 0.5*w**2*r_squared*(1-self.alpha**2) + 2*self.alpha*w
        
        #Coulomb repulsio
        r_12 = system.particle_distance_squared()
        #energy += 1/r_12
            #old_position = particle.position
            #particle.random_move(old_position)
                
        return energy
    def local_energy2(self, r, wf, system):
        #local energy with coulomb interaction
        total_energy = 0.0
        #2 electron case, ground state, opposite spins, a= 1.0
        a = 1.0
        for i in range(system.number_of_particles):
            r_single_particle = 0.0
            for j in range(system.dimensions):
                r_single_particle += r[i,j]**2 #Why not just r[i,j]?
            
            r_12 = math.sqrt(r_single_particle)
            denom = a/((1+self.beta*r_12)**2)
            total_energy += 2*self.alpha**2*system.w*r_single_particle
            total_energy -= 4*self.alpha*system.w
            total_energy -= 2*self.alpha*system.w*denom
            total_energy += 2*denom*(a*denom + 1/r_12 - 2*self.beta/(1+self.beta*r_12))
        return total_energy
        
    def relative_distance(self, system,r):
        n_particles = system.number_of_particles
        dimensions = system.dimensions
        for i1 in range(n_particles-1):
            for i2 in range(i1+1, n_particles):
                r_12 = 0.0
                for j in range(dimensions):
                    r_12 += (r[i1,j] - r[i2,j])**2
                    
                return r_12
        
    def MC_calculations(self,system, energy_min = 100.0):
        n_particles = system.number_of_particles
        dimensions = system.dimensions
        step_length = self.step_length
        n_cycles = self.mc_cycles
        #variations = self.variations
#        r_0 = np.zeros((n_particles,dimensions), np.double)
#        r_n = np.zeros((n_particles,dimensions), np.double)
        
        energy = energy2 = 0.0
        accept = 0.0
        delta_e = 0.0
        #Initial position
#        for i in range(n_particles):
#            for j in range(dimensions):
#                r_0[i,j] = step_length * (random() - .5)
        #Initial position
        system.advance_time(step_length)
        wf_0 = self.trial_wavefunction(r_0, system)
        for cycle in range(self.mc_cycles):
            #Trial position
            system.advance_time(step_length)
#            for i in range(n_particles):
#                for j in range(dimensions):
#                    r_n[i,j] = r_0[i,j] + step_length * (random() - .5)
                    #r_n[i,j] = r_0[i,j] + step_length * (random.normal_distribution(0.0,1.0) - .5)
    
            wf_n = self.trial_wavefunction(r_n, system)
            
            #Metropolis test to see whether we accept the move
            if random() < wf_n**2 / wf_0**2:
                r_0 = r_n.copy()
                wf_0 = wf_n
                accept += 1
            #update expectation values
            
            delta_e = self.local_energy(r_0, wf_0, system)
            #delta_e = self.local_energy2(r_0, wf_0, system)
            energy += delta_e
            energy2 += delta_e**2
        #We calculate mean, variance and error ...
        energy /= n_cycles
        energy2 /= n_cycles
        variance = abs(energy2 - energy**2)
        error = math.sqrt(variance/n_cycles);
        results = self.alpha, self.beta, energy, variance, error, accept*1.0/n_cycles
        return results
        
            
            
    def write_to_file(self,outfile, results):    
        #print(results)
        outfile.write('%f %f %f %f %f %f\n' %(results[0],results[1],results[2],results[3],results[4],results[5]))
            
            #outfile.write('%f %f %f %f %f %f\n' %(alpha,beta,avg_energy,variance,error,accept*1.0/(n_cycles)))
            #...and write them to file
#        avg_energy = sum_energy/variations
#        if avg_energy < energy_min:
#            energy_min = avg_energy
#            alpha_min = solver.alpha
#            beta_min = beta            

#for alpha in alphas:
#    for beta in betas:
#        sum_energy = 0.0
#        print("Calculating alpha = {}, beta = {}...".format(alpha,beta))
#        for variation in range(max_variations):
#            #alpha += alpha_step
#        
#            energy = energy2 = 0.0
#            accept = 0.0
#            delta_e = 0.0
#        
#            
#            #Loop over MC cycles

#        #
            #We calculate mean, variance and error ...
#            energy /= n_cycles
#            sum_energy += energy
#            energy2 /= n_cycles
#            variance = abs(energy2 - energy**2)
#            #print(energy, energy2, variance)
#            #print(variance,type(variance))
#            error = math.sqrt(variance/n_cycles);
#            #...and write them to file
#        avg_energy = sum_energy/max_variations
#        if avg_energy < energy_min:
#            energy_min = avg_energy
#            alpha_min = alpha
#            beta_min = beta
#        outfile.write('%f %f %f %f %f %f\n' %(alpha,beta,avg_energy,variance,error,accept*1.0/(n_cycles)))
#        
#outfile.close()
#print("Minimum energy {} at alpha = {}, beta = {}".format(energy_min,alpha_min, beta_min))
#print('\nDone. Results are in the file "%s", formatted as:\n\
#alpha, beta, <energy>, variance, error, acceptance ratio' %(outfilename))
#t_f = time.time()
#delta_t = t_f - t_i
#print("computation time took {}s".format(delta_t))
#
#            