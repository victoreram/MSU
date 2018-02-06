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
    def __init__(self, dimensions=2, step_length=1.0):
        #dimensions = system.dimensions
        self.step_length = step_length
        self.position = np.zeros(dimensions)
        
    def random_move(self, step_length=1.0):
        '''Brute Force move; randomly choose where to move without importance
        sampling'''
        #print(self.position)
        temp_coordinates = self.position.copy()
        for index in range(len(self.position)):
            old_coordinate = self.position[index]
            new_coordinate = old_coordinate + step_length * (random() - .5)
            temp_coordinates[index] = new_coordinate
            #self.position[index] = new_coordinate
        return temp_coordinates
        #print(self.position)
    def move(self, step_length=1.0):
        '''Move for importance sampling'''
        pass
    
    def distance_difference(self, old_position):
        arg = 0.0
        for i in range(len(old_position)):
            change = self.position[i] - old_position[i]
            arg += change**2
            
        total = math.sqrt(arg)
        return total
        
        
    def r_squared(self):
        r_sum = 0.0
        for r in self.position:
            r_sum += r**2
        return r_sum