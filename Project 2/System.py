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
        #self.particles = [Particle(self) for particle in range(dimensions)]
        
    def __str__(self):
        string = '''***SYSTEM PARAMETERS***\nNumber of particles: {}\nDimensions: {}\nOscillator frequency w: {}
        '''.format(self.number_of_particles, self.dimensions, self.w)
        return string
