class Planet():
    def __init__(self, name="Planet X", mass=1, coordinates=[1,0,0], velocity=[0,0,0]):
        self.mass = mass
        self.name = name
        self.coordinates = coordinates
        self.x = coordinates[0]
        self.y = coordinates[1]
        self.z = coordinates[2]
        self.velocity = velocity
        self.vx = velocity[0]
        self.vy = velocity[1]
        self.vz = velocity[2]
        self.r = self.orbit_radius()
        self.velocity_mag = (self.vx**2 + self.vy**2 + self.vz**2)**(0.5)
        self.position = (self.x,self.y,self.z)

    def orbit_radius(self):
        return (self.x**2+self.y**2+self.z**2)**(0.5)
    
    def gravitational_acceleration(self, other_planet, coord_index):
        difference = self.coordinates[coord_index]-other_planet.coordinates[coord_index]
        return -4*m.pi*self.mass*difference/((self.distance(other_planet))**3)
    
    def distance(self, other_planet):
        return ((self.x-other_planet.x)**2 + (self.y-other_planet.y)**2 + (self.z-other_planet.z)**2)**(0.5)
    
    def kinetic_energy(self):
        return (.5*self.mass*((self.vx)**2 + (self.vy)**2 + (selfv.z)**2))
    
    def potential_energy(self, other_planet):
        G=1
        return G*self.mass*other_planet.mass/(self.distance(other_planet))
        
    def angular_momentum(self, other_planet):
        return mass*self.distance(other_planet)*self.velocity_mag
    
    