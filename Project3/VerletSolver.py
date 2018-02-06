import math
class VerletSolver():
    def __init__(self, h, time):
        
        self.h = h
        self.time = time
        
    def coordinate(xi,h,vi,ai):
        return xi + h*vi + h**2*ai/2

    def velocity(vi,h,a_i_1, a_i):
        return vi + (h/2)*(a_i_1+a_i)

    def acceleration(coord, dist):
        if dist == 0:
            return 0
        else:
            return -4*coord*math.pi**2/(dist**3)
    
    def kinetic_energy(m, v):
         return (.5*m*v**2)
    
    def potential_energy(m1, m2, d, G=6.67E-11):
        return G*m1*m2/d
    
    def angular_momentum(m, v, d):
        return m*v*d
    
    def algo(self, c_i, v_i, r, h):
        a_i = self.acceleration(c_i, r)
        c_i_1 = self.coordinate(c_i, h, v_i, a_i)
        a_i_1 = self.acceleration(c_i_1, r)
        v_i_1 = self.velocity(v_i, a_i_1, a_i)
        return c_i_1, v_i_1
        
        
        