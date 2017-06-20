# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 14:02:04 2017

@author: ramir
"""
import numpy as np
from random import random
import math
#from System import System

class Particle:
    def __init__(self, dimensions):
        #dimensions = system.dimensions
        self.position = np.zeros(dimensions)
        
    def random_move(self, step_length=1.0):
        '''Brute Force move; randomly choose where to move without importance
        sampling'''
        #print(self.position)
        for index in range(len(self.position)):
            old_coordinate = self.position[index]
            new_coordinate = old_coordinate + step_length * (random() - .5)
            self.position[index] = new_coordinate
            
        #print(self.position)
    def move(self, step_length):
        
        pass
    
    def distance_difference(self, old_position):
        arg = 0.0
        for i in range(len(old_position)):
            change = self.position[i] - old_position[i]
            arg += change**2
            
        total = math.sqrt(arg)
        return total
        
    def relative_distance_squared(self):
        r_squared = 0.0
        for coordinate in self.position:
            r_squared += coordinate**2
        return r_squared
        
    def r_squared(self):
        r_sum = 0.0
        for r in self.position:
            r_sum += r**2
        return r_sum

            
            
#sys = System(2,2,1.0)
e = Particle(2)
e.random_move(1.0)