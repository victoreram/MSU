# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:42:20 2017

@author: ramir
"""
import math
import numpy as np
from System import System
from Particle import Particle
from random import random
class Solver:
    def __init__(self, system, mc_cycles, alpha, beta, alpha_step,beta_step,alpha_variations,beta_variations):
        '''Solver(mc_cycles, alpha, beta, alpha step, beta step, 
        number of alpha variations, numbe of beta variations'''
        self.system = system        
        self.alpha_variations = alpha_variations
        self.beta_variations = beta_variations
        self.mc_cycles = mc_cycles
        self.alpha = alpha
        self.beta = beta
        #self.hamiltonian = hamiltonian
        #range of parameters to check; doesn't change with changing alpha
        self.alpha_step = alpha_step
        self.beta_step = beta_step
        self.final_alpha = alpha + (alpha_variations-1)*alpha_step
        self.final_beta = beta + (beta_variations-1)*beta_step
        self.alphas = np.linspace(alpha,self.final_alpha,alpha_variations)
        self.betas = np.linspace(beta,self.final_beta,beta_variations)
    #precalculated parameters that update with changing alpha, omega
    @property 
    def two_alpha_w(self):
        return 2*self.alpha*self.system.w
    
    @property
    def one_minus_alpha_squared(self):
        return 1-self.alpha**2
        
    @property
    def half_omega_squared(self):
        return 0.5*self.system.w**2
        
    def optimize_parameters(self, outfile, minimum_energy=1E9):
        '''Test varational parameters. Default parameters include alpha
        and beta'''
        alphas = self.alphas.copy()
        betas = self.betas.copy()
        number_of_particles = self.system.number_of_particles
        dimensions = self.system.dimensions
        w = self.system.w
        step_length = self.system.step_length
        for alpha in alphas:
            for beta in betas:
                #change variational parameters of solver
                self.alpha = alpha
                self.beta = beta
                #reinitialize particles in system to beginning postion after each set of alpha and beta               
                self.system = System(number_of_particles,dimensions,w,step_length)
#                self.system.particles = [Particle(dimensions, step_length) for p in range(number_of_particles)]
#                self.system.position_matrix = [np.zeros(dimensions) for particle in range(number_of_particles)]
                #print("Inintial relative distance = ", self.system.position_matrix)
                print("calculating alpha = {}, beta = {}...".format(alpha,beta))
                results = self.MC_calculations()
                self.write_to_file(outfile,results)
                energy = results[2]
                #find mininmum
                if energy < minimum_energy:
                    minimum_energy = energy
                    minimum_alpha = results[0]
                    minimum_beta = results[1]
        print("Minimum energy {} at alpha = {}, beta = {}".format(minimum_energy, minimum_alpha, minimum_beta))
        return minimum_alpha, minimum_beta, minimum_energy                      

    def trial_wavefunction(self, r):
        wf = 1.0
        r_sum = 0.0
        for i in range(self.system.number_of_particles):
            #loop over each particle
            r_ij_particle = 0.0
            for j in range(self.system.dimensions):
                r_ij_particle += r[i][j]**2
            r_sum += r_ij_particle
        wf *= math.exp(-0.5*self.alpha*self.system.w*r_sum)
#        for p in range(self.system.number_of_particles):
#            r_sum = self.system.relative_distance_squared(p)
#            wf *= math.exp(-0.5*self.alpha*self.system.w*r_sum)
        #Jastrow Factor
        #wf = self.jastrow(wf)
        #arg = self.system.particle_distance_squared(r_sum)
        #wf *= math.exp(0.5*arg/(1.0+self.beta*arg))
        return wf
        
    def jastrow(self,wf):
        for i in range(self.system.number_of_particles-1):
            for j in range(i+1,self.system.number_of_particles):
                r_12 = 0.0
                for k in range(self.system.dimensions):
                    r_12 += (self.system.position_matrix[i][k] - self.system.position_matrix[j][k])**2
                arg = math.sqrt(r_12)
                wf *= math.exp(0.5*arg/(1.0+self.beta*arg))
        return wf
 
        
    def local_energy(self, r, a = 1.0):
        '''Local Energy using analytical expression. '''
        w = self.system.w
        alpha = self.alpha
        beta = self.beta
        energy = 0.0
        #energy += self.two_alpha_w
        #r_1_squared_plus_r_2_squared = self.system.r_sum_squared()
        r_sum = 0.0
        for i in range(self.system.number_of_particles):
            #loop over each particle
            r_ij_particle = 0.0
            for j in range(self.system.dimensions):
                r_ij_particle += r[i][j]**2
            r_sum += r_ij_particle
        #print("Particle coordinates: {}, r_1^2 + r_2^2 = {} ".format( self.system.position_matrix, r_1_squared_plus_r_2_squared))
        #r_12 = math.sqrt(self.system.particle_distance_squared())
        #r_12_inverse = 1/r_12
        #one_plus_beta_r_12 = 1+self.beta*r_12
        #calculation of energy from A.1.18 of Christian's Thesis
        energy += self.half_omega_squared*r_sum*self.one_minus_alpha_squared
        #print("delta e = {}*{}*{} = {}".format(self.half_omega_squared,r_1_squared_plus_r_2_squared,self.one_minus_alpha_squared,energy))
        #energy -= 2*a*beta/(one_plus_beta_r_12**3)
        #energy -= a**2/(one_plus_beta_r_12**4)
        #energy += (alpha*w*r_12 + r_12_inverse)*a/(one_plus_beta_r_12**2)
        #Coulomb repulsion        
        #energy += r_12_inverse
        #print("Energy = 0.5*{}^2*{}*(1-{}^2) + 2*{}*w = {}".format(w,round(r_squared,2),self.alpha,self.alpha,energy))
        #r_12 = system.particle_distance_squared()
        #energy += 1/r_12
        return energy
        

    def MC_calculations(self, energy_min = 1.0E9):
        n_cycles = self.mc_cycles
        #r_0 = np.zeros((n_particles,dimensions), np.double)
        #r_n = np.zeros((n_particles,dimensions), np.double)
        energy = 0.0
        energy2 = 0.0
        accept = 0.0
        delta_e = 0.0
        #iters=0
#        for i in range(n_particles):
#            for j in range(dimensions):
#                r_0[i,j] = step_length * (random() - .5)
        #Initial position
        #print("Initial position: ", self.system.position_matrix)

        r_0 = self.system.advance_time()
        wf_0 = self.trial_wavefunction(r_0)
        #energy += self.two_alpha_w*self.mc_cycles
        #print("Initial wavefunction = {}".format(wf_0))
        for cycle in range(self.mc_cycles):
            #Trial position
            r_n = self.system.advance_time()
            #print("Particle coordinates: ", self.system.position_matrix)
            wf_n = self.trial_wavefunction(r_n)
            #Metropolis test to see whether we accept the move
            #get ZeroDivisionError without "if wf_0 == 0.0" <- resolved as of 9/14
            roll = random()            
            #print("random: {:.3} new wf: {:.3} old wf: {:.3} ratio: {:.3}".format(roll, wf_n**2, wf_0**2, wf_n**2 / wf_0**2))
            if roll < wf_n**2 / wf_0**2:
                wf_0 = wf_n
                accept += 1
                self.system.position_matrix = r_n
                r_0 = r_n
                #print("Move accepted; moving to ", r_n)

                
            #update expectation values            
            delta_e = self.local_energy(r_0) + self.two_alpha_w
            energy += delta_e
            energy2 += delta_e**2
            #print("{} + {} = {}".format(old_energy,delta_e,energy))

            #iters += 1
        #We calculate mean, variance and error ...
        energy /= n_cycles
        energy2 /= n_cycles
        print("Energy = ", energy)
        variance = abs(energy2 - energy**2)
        error = math.sqrt(variance/n_cycles);
        results = self.alpha, self.beta, energy, variance, error, accept*1.0/n_cycles
        print("variance: {:e} error: {:e} acceptance ratio: {}".format( variance, error, accept*1.0/n_cycles)) 
        return results
            
    def write_to_file(self,outfile,results):    
        outfile.write('%f %f %f %f %f %f\n' %(results[0],results[1],results[2],results[3],results[4],results[5]))
            
    def quantum_force(self):
        
#        {
#  // compute the first derivative 
#  for (int i = 0; i < NumberParticles; i++) {
#    for (int k = 0; k < Dimension; k++) {
#      // single-particle part, replace with Slater det for larger systems
#      double sppart = DerivativeSPwavefunction(r(i,k),VariationalParameters(0));
#      //  Jastrow factor contribution
#      double Jsum = 0.0;
#      for (int j = 0; j < NumberParticles; j++) {
#	if ( j != i) {
#	  Jsum += JastrowDerivative(r, VariationalParameters(1), i, j, k);
#	}
#      }
#      qforce(i,k) = 2.0*(Jsum+sppart);
#    }
#  }
#} // end of QuantumForce function
        
        '''
        Calculate 1st derivative
        Apply Jastrow Factor Contribution using Jastrow Derivative
        '''
        pass
    
    def Jastrow_derivative(self):
        pass

#inline double JastrowDerivative(Matrix &r, double beta, int i, int j, int k){
#  return (r(i,k)-r(j,k))/(RelativeDistance(r, i, j)*pow(1.0+beta*RelativeDistance(r, i, j),2));

        
#    def local_energy2(self, r, wf, system):
#        #local energy with coulomb interaction
#        total_energy = 0.0
#        #2 electron case, ground state, opposite spins, a= 1.0
#        a = 1.0
#        for i in range(system.number_of_particles):
#            r_single_particle = 0.0
#            for j in range(system.dimensions):
#                r_single_particle += r[i,j]**2 #Why not just r[i,j]?
#            
#            r_12 = math.sqrt(r_single_particle)
#            denom = a/((1+self.beta*r_12)**2)
#            total_energy += 2*self.alpha**2*system.w*r_single_particle
#            total_energy -= 4*self.alpha*system.w
#            total_energy -= 2*self.alpha*system.w*denom
#            total_energy += 2*denom*(a*denom + 1/r_12 - 2*self.beta/(1+self.beta*r_12))
#        return total_energy
#        
#    def relative_distance(self, r):
#        n_particles = self.system.number_of_particles
#        dimensions = self.system.dimensions
#        for i1 in range(n_particles-1):
#            for i2 in range(i1+1, n_particles):
#                r_12 = 0.0
#                for j in range(dimensions):
#                    r_12 += (r[i1,j] - r[i2,j])**2
#                    
#                return r_12
        
            
###############################################################################            
#    def trial_wavefunction(self, r, system):
#        r_sum = 0.0
#        for i in range(system.number_of_particles):
#            #loop over each particle
#            r_ij_particle = 0.0
#            for j in range(system.dimensions):
#                r_ij_particle += r[i,j]**2
#            r_sum += r_ij_particle
#        wf = math.exp(-0.5*self.alpha*system.w*r_sum)
#        if self.jastrow_bool == True:
#            wf = self.jastrow(wf, r, system)
#        return wf
#        
        
#    def jastrow(self,wf,r,system):
#        for i in range(system.number_of_particles-1):
#            for j in range(i+1,system.number_of_particles):
#                r_12 = 0.0
#                for k in range(system.dimensions):
#                    r_12 += (r[i,k] - r[j,k])**2
#                arg = math.sqrt(r_12)
#                wf *= math.exp(0.5*arg/(1.0+self.beta*arg))
#        return wf
        #analytic expression version
        #for i in range(system.number_of_particles)
        
#    def local_energy(self,r, wf, system, h=0.001, h2 = 1E6):
#        '''Local energy using numerical derivative'''
#        #Kinetic energy
#        n_particles = system.number_of_particles
#        dimensions = system.dimensions
#        r_plus = r.copy()
#        r_minus = r.copy()
#        e_kinetic = 0.0
#        #numerical derivative
#        for i in range(n_particles):
#            for j in range(dimensions):
#                r_plus[i,j] = r[i,j] + h
#                r_minus[i,j] = r[i,j] - h
#                wf_minus = self.trial_wavefunction(r_minus, system)
#                wf_plus = self.trial_wavefunction(r_plus, system)
#                e_kinetic -= wf_minus+wf_plus-2*wf;
#                r_plus[i,j] = r[i,j]
#                r_minus[i,j] = r[i,j]
#        
#        e_kinetic = .5*h2*e_kinetic/wf
#        #Potential energy
#        e_potential = 0.0
#        #harmonic oscillator  contribution
#        for i in range(n_particles):
#            r_single_particle = 0.0
#            for j in range(dimensions):
#                r_single_particle += r[i,j]**2
#            e_potential += 0.5*r_single_particle
#    
#        #Electron-electron contribution
#        for i1 in range(n_particles-1):
#            for i2 in range(i1+1,n_particles):
#                r_12 = 0.0
#                for j in range(dimensions):
#                    r_12 += (r[i1,j] - r[i2,j])**2
#                if self.jastrow_bool == True:
#                    e_potential += 1.0/math.sqrt(r_12)
#        
#        return e_potential + e_kinetic


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