class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def __mul__(self, v2):
        x2 = self.x*v2.x
        y2 = self.y*v2.y
        return Vector(x2, y2)
    
    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"
    
    def __add__(self, v2):
        x2 = self.x + v2.x
        y2 = self.x + v2.y
        return Vector(x2,y2)
    
    def __sub__(self, v2):
        x2 = self.x - v2.x
        y2 = self.x - v2.y
        return Vector(x2,y2)
    
    def dot(self, v2):
        return self.x*v2.x + self.y*v2.y
        
    def magnitude(self):
        return (self.x**2+self.y**2)**(0.5)