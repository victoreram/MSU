# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:39:40 2017

@author: ramir
"""
from Particle import Particle
class System:
    def __init__(self,number_of_particles,dimensions,w):
        '''System(number_of_particles, dimensions, oscillator frequency w))'''
        self.number_of_particles = number_of_particles
        self.dimensions = dimensions
        self.w = w
        #Generate particles and store them in a list
        self.particles = [Particle(self.dimensions) for p in range(number_of_particles)]
        #Store each particle's position in a matrix        
        self.position_matrix = [particle.position for particle in self.particles]
        
    def __str__(self):
        string = '''***SYSTEM PARAMETERS***\nNumber of particles: {}\nDimensions: {}\nOscillator frequency w: {}
        '''.format(self.number_of_particles, self.dimensions, self.w)
        return string
        
    def particle_distance_squared(self):
        n_particles = self.number_of_particles
        r_12 = 0.0
        for i1 in range(n_particles-1):
            for i2 in range(i1+1,n_particles):
                #r_12 = 0.0 here originally
                particle_1 = self.particles[i1]
                particle_2  = self.particles[i2]
                r_12 += particle_1.distance_difference(particle_2.position)
        return r_12
        
    def relative_distance_squared(self):
        
        r_sum = 0.0
        for i in range(self.number_of_particles):
            #loop over each particle
            r_ij_particle = 0.0
            for j in range(self.dimensions):
                r_ij_particle += self.position_matrix[i][j]**2
            r_sum += r_ij_particle
            
        return r_sum
        
    def advance_time(self,step_length):
        for particle in self.particles:
            particle.move(step_length)
        