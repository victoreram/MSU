# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 14:02:04 2017

@author: ramir
"""
import numpy as np
from random import random
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
            
#sys = System(2,2,1.0)
e = Particle(2)
e.random_move(1.0)