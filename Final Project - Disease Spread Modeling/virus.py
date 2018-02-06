class Virus():
    
    def __init__(self, x, y, grid, infectability = 0.1, lethality = 0.1, mutationrate = 0.1, immunity = 0.5, inhost = False):
        '''Constructs instance of virus with the following attributes:
        infectability: how infectuous the virus is
        lethality: the virus's ability to kill
        mutationrate: chance that virus mutates
        lifetime: how many iterations does the virus live
        inhost: if True, this virus is affecting a person'''
        self.x = x
        self.y = y
        self.grid = grid
        self.infectability = infectability
        self.lethality = lethality
        self.immunity = immunity
        self.mutationrate = mutationrate
        self.inhost = inhost
        
    def check_surroundings(self):
        '''Method used to check surroundings. '''
        surroundings = {}

        for x in range(self.x - 1, self.x + 2):
            for y in range(self.y - 1, self.y + 2):
                if (not (x == self.x and y == self.y)) and self.on_grid(x, y):
                    
                    # 1 is a placeholder for person; ideally I could access an instance of a person in this position
                    if self.grid(x,y) == 1:
                        near_person = True
                        person_x = x
                        person_y = y
                        
                        
        return [person_x, person_y, near_person, person_instance]
        
        
    def on_grid(self, x, y):
        '''Checks if coordinates given (x,y) are valid locations of the grid's dimensions (r,c)'''
        return x >= 0 and y >= 0 and x < self.grid.c and y < self.grid.r
    
    def infect(self):
        '''Checks surroundings for a person to infect. If a person is found, roll to see if infection is successful.
        If so, change instance of that person to be "infected"'''
        result = self.check_surroundings()
        if result[2]:
            roll = random.uniform(0,1)
            if roll < self.infectability:
                # infect the instance of a person in this position
                coordinates = (result[0],result[1])
                # p = person_dict[coordinates] <- access the person
                # p.infected = True <- change state of person
            
    def move(self):
        '''looks for open position in grid and assigns virus's x and y coordinates to that position'''
        open_pos = []
        for x in range(self.x - 1, self.x + 2):
            for y in range(self.y - 1, self.y + 2):
                if not (x == self.x and y == self.y):
                    # Allow the virus to divide to the cell at x, y if the point
                    # is on the grid and no existing virus lives at x, y.
                    if self.on_grid(x, y) and self.grid.bac_mat[y][x] == 0:
                        open_pos.append([x,y])
                        
        random.shuffle(open_pos)
        #print(open_pos)
        self.x = open_pos[0][0]
        self.y = open_pos[0][1]


    
    def mutate(self):
        '''Each immunity gene has a [rate] chance of mutating'''
        roll = random.uniform(0,1)
        if self.mutationrate > roll:
            new_attributes = [random.uniform(0,1) for _ in range(4)]
            self.infectability = new_attributes[0]
            self.lethality = new_attributes[1]
            self.mutationrate = new_attributes[2]
            self.immunity = new_attributes[3]
        #self.imm = [random.uniform(0, 1) if random.random() < self.mutationrate else current_val for current_val in self.imm]
        return self
        
    def die(self):
        '''Check if virus is immune to anti-virus; if virus doesn't pass check, it dies.'''   
        for i,j in zip(self.immunity, self.grid.get_anti_bac(self.x)):
            if j > i:
                self.grid.bac_mat[self.y][self.x] = 0
                del self.grid.bac[(self.y, self.x)]
                return True
                
        return False
    
    
    def infected_iterate(self):
        '''Method that describes a virus actions inside the body. 
        Approach in one of two ways:
        1) virus stays "inside" person; i.e. shares location of the person, checks surroundings and infects surrounding
        2) delete instance of this virus, and all of the person-infecting-other-people interactions are within the 
        person class
        '''
        
        #1)
        #self.x = person.x
        #self.y = person.y
        result = check.surroundings()
        self.infect()
        
        #in case we account for dying people
        roll = random.uniform(0,1)
        if roll > self.lethality:
            self.grid.bac_mat[self.y][self.x] = 0
            del self.grid.bac[(self.y, self.x)] #placeholder for deleting person
        
        self.mutate()
        
    
    def iterate(self):
        '''
        How virus interacts;
        If the virus is in the host, it acts thru the "infected_iterate" method
        If the virus is not, it keeps moving and finds a person to infect
        '''
        
        if self.inhost == True:
            self.infected_iterate()
            
        else:
            self.move()
            self.infect()
            self.mutate()
        


    def draw(self):
        '''Draws the virus and represents the colors based on the immunities'''
        #rgb = self.imm
        plt.scatter(self.y, self.x)    
    
    