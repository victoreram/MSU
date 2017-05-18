import math

class Vector():
    def __init__(self,x=0,y=0,z=0):
        self.x = x
        self.y = y
        self.z = z
        self.components = [self.x,self.y,self.z]
        
    def __add__(self, v2):
        return Vector(self.x+v2.x, self.y+v2.y, self.z+v2.z)
    
    def __sub__(self, v2):
        return Vector(self.x-v2.x, self.y-v2.y, self.z-v2.z)
    
    def __mul__(self, n):
        if isinstance(n, Vector):
            return Vector(self.x*n.x, self.y*n.y, self.z*n.z)
        elif isinstance(n, int) or isinstance(n, float):
            return Vector(self.x*n, self.y*n, self.z*n)
        
    def __div__(self,n):
        if isinstance(n, Vector):
            return Vector(self.x/n.x, self.y/n.y, self.z/n.z)
        elif isinstance(n, int) or isinstance(n, float):
            return Vector(self.x/n, self.y/n, self.z/n)
        
        
    def __str__(self):
        return "({},{},{})".format(self.x,self.y,self.z)
    
    def distance(self, v2):
        return ((self.x-v2.x)**2 + (self.y-v2.y)**2 + (self.z-v2.z)**2)**0.5
    
    def dot(self, v2):
        return self.x*v2.x + self.y*v2.y + self.z*v2.z
    
    def magnitude(self):

        mag_sq = 0
        for comp in self.components:
            mag_sq += comp**2
        return math.sqrt(mag_sq)
    
    def distance(self, v2):
        return ((self.x-v2.x)**2 + (self.y-v2.y)**2 +(self.z-v2.z)**2)**(0.5)
        