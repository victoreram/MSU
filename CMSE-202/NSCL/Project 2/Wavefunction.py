# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:01:59 2017

@author: ramir
"""

class Wavefunction:
    def __init__(self, parameters, system):
        self.system = system
        self.parameters = parameters
        self.wavefunctions = {[False, False] : self.two_electron()}
        self.value = self.wavefunctions[parameters]
        
    def two_electron(self):
        '''Default:
        -Two Electron
        -Importance Sampling Disabled
        -Jastrow Factor Disabled'''
        pass
    
    def two_electron_with_importance(self):
        pass
    
    def two_electron_random_with_jastrow(self):
        pass
    
    def two_electron_with_jastrow_coulomb(self):
        pass
