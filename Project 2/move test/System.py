# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:39:40 2017

@author: ramir
"""
from Particle import Particle
class System:
    def __init__(self,number_of_particles,dimensions,w,step_length):
        '''System(number_of_particles, dimensions, oscillator frequency w, step length))'''
        self.number_of_particles = number_of_particles
        self.dimensions = dimensions
        self.w = w
        self.step_length = step_length
        #Generate particles and store them in a list
        self.particles = [Particle(self.dimensions, self.step_length) for p in range(self.number_of_particles)]
        #Store each particle's position in a matrix 
        self.position_matrix = [particle.position for particle in self.particles]
        self.greens_function = 1.0
#    for automatically updating attributes
#    @property
#    def position_matrix(self):
#        return [particle.position for particle in self.particles]
        
    def __str__(self):
        string = '''***SYSTEM PARAMETERS***\nNumber of particles: {}\nDimensions: {}\nOscillator frequency w: {}
        '''.format(self.number_of_particles, self.dimensions, self.w)
        return string
        
    def particle_distance_squared(self, r_12 = 0.0):
        '''Returns sum of distances between particles squared, aka r_12^2'''
        n_particles = self.number_of_particles
        for i1 in range(n_particles-1):
            for i2 in range(i1+1,n_particles):
                #r_12 = 0.0 here originally
                particle_1 = self.particles[i1]
                particle_2  = self.particles[i2]
                for j in range(self.dimensions):
                    r_12 += (particle_1.position[j]-particle_2.position[j])**2

                #r_12 += particle_1.distance_difference(particle_2.position)
        return r_12
        
    def r_sum_squared(self):
        '''Returns sum of relative distances among all particles squared, aka r_1^2 + r_2^2 + ...'''
        r_sum = 0.0
        for i in range(self.number_of_particles):
            #loop over each particle
            r_ij_particle = 0.0
            for j in range(self.dimensions):
                r_ij_particle += self.position_matrix[i][j]**2
            r_sum += r_ij_particle
#        for j in range(self.dimensions):
#            r_sum += self.position_matrix[i][j]**2
            
        return r_sum
        
    def advance_time(self):
        '''Moves each particle in system within a defined range by step length'''
#        for particle in self.particles:
#            particle.random_move(self.step_length)
        '''Move one particle at a time'''
        temp_position_matrix = self.position_matrix.copy()
        for p in range(len(self.particles)):
            temp_position_matrix[p] = self.particles[p].random_move(self.step_length)
            
        return temp_position_matrix
        #self.position_matrix = [particle.position for particle in self.particles]
        #print(self.position_matrix)

        