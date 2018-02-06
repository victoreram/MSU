# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 22:55:17 2017

@author: ramir
"""

class test:
    def __init__(self, a, b):
        self.a = a
        self.b = b
    @property
    def c(self):
        return self.a + self.b        
test1 = test(2,3)
print(test1.c)
test1.a = 3
print(test1.c)