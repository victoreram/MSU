# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:01:59 2017

@author: ramir
"""
#Hamiltonians = {[False, False] : two_electron(system)}
class Hamiltonian:
    def __init__(self, parameters, system):
        '''Hamiltonian class used to calculate the Hamiltonian given a set
        of parameters. 
        (Type) Attribute:
            (System) system: Contains parameters such as harmonic oscillator potential,
            number of particles, and dimensions
            (List) parameters: Dictate which factors are enabled as a list of booleans
            Structure is as follows:
                [jastrow_factor_enabled, importance_sampling_enabled]
            (Dictionary) hamiltonians: Contain different formulae to calculate Hamiltonians
                Key:Value = parameters:hamiltonian'''
            
        self.system = system
        self.parameters = parameters
        self.hamiltonians = {[False, False] : self.two_electron(),
                             [True, False] : self.two_electron_random_with_jastrow()}
        #self.value = self.Hamiltonians[parameters]
        
    def two_electron(self):
        '''Default:
        -Two Electron
        -Importance Sampling Disabled
        -Jastrow Factor Disabled'''
        '''Local Energy using analytical expression. '''
        w = self.system.w
        energy = 0.0
        r_squared = 0.0
        for particle in self.system.particles:
            r_i = particle.r_squared()
            r_squared += r_i
        energy += 0.5*w**2*r_squared*(1-self.system.alpha**2) + 2*self.system.alpha*w
        #print("Energy = 0.5*{}^2*{}*(1-{}^2) + 2*{}*w = {}".format(w,round(r_squared,2),self.alpha,self.alpha,energy))
        return energy
    
    def two_electron_with_importance(self):
        pass
    
    def two_electron_random_with_jastrow(self):
        '''
        -Two Electron
        -Importance Sampling Disabled
        -Jastrow Factor Disabled'''
        '''Local Energy using analytical expression. '''
        w = self.system.w
        energy = 0.0
        r_squared = 0.0
        for particle in self.system.particles:
            r_i = particle.r_squared()
            r_squared += r_i
        energy += 0.5*w**2*r_squared*(1-self.system.alpha**2) + 2*self.system.alpha*w
        #print("Energy = 0.5*{}^2*{}*(1-{}^2) + 2*{}*w = {}".format(w,round(r_squared,2),self.alpha,self.alpha,energy))
        #Coulomb repulsion
        r_12 = self.system.particle_distance_squared()
        energy += 1/r_12
        #old_position = particle.position
        #particle.random_move(old_position)
        return energy
        pass

def two_electron_with_jastrow_coulomb(self):
    pass
