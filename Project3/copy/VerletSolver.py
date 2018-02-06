import math as m

class VerletSolver():
    def __init__(self, h, time):
        
        self.h = h
        self.time = time
        
    def coordinate(xi,h,vi,ai):
        return xi + h*vi + h**2*ai/2

    def velocity(vi,h,a_i_1, a_i):
        return vi + (h/2)*(a_i_1+a_i)

    def acceleration(coord, dist):
        return -4*coord*m.pi**2/(dist**3)
    
    def kinetic_energy(planet):
        return (.5*planet.mass*((planet.vx)**2 + (planet.vy)**2 + (planet.vz)**2))
    
    def potential_energy(planet, other_planet):
        G=1
        return G*planet.mass*other_planet.mass/(planet.distance(other_planet))
    
    def angular_momentum(planet, other_planet):
        return mass*planet.distance(other_planet)*planet.velocity_mag
        
        