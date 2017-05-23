import numpy as np
from sympy import hermite
from sympy.abc import x
import scipy
from Vector import Vector as vec
import math
import scipy.linalg as la
from decimal import *
#import pandas as pd


def magic(n):
    #given shell number n, returns magic number at n
    result = 0
    for i in range(1,n+1):
        result += 2*i
    return result

def energy(n):
    #given shell number n, returns energy at n
    result = 0
    for i in range(1,n+1):
        result += 2*i**2
    return result

def nx_ny(n):
    #given shell number n, returns tuple of possible quantum numbers in (nx, ny) form
    quantum_numbers = []
    for i in range(n+1):
        quantum_numbers.append(vec(i,n-i))
    return tuple(quantum_numbers)

def n_m(e):
    quantum_numbers = []
    for m in range(e-1,-e,-2):
        quantum_numbers.append(vec((e-abs(m)-1)//2,m))
    #for m in range(-e,e,2):
    #    quantum_numbers.append((()))
    #for m in range(e-1,-1,-2):
        #quantum_numbers.append(((e-abs(m)-1)//2,m))
        #quantum_numbers.append(((e-abs(m)-1)//2,-m))


    return tuple(quantum_numbers)

def wavefunction(n):
    #given shell number n, returns the hermite polynomial at that n
    return hermite(n,x)


class single_particle():

    def __init__(self, i, n, m, m_s, state):
        #takes in n, corresponding to shell number, where n = 1,2,3...
        self.shell = i
        self.n = n
        self.m = m
        self.degeneracy = 2*i
        self.E_nx_ny = str(i)
        self.magic_number = magic(i)
        self.energy = energy(i)
        self.nx_ny = []
        self.n_m = [n,m]
        self.state = state
        self.m_s = m_s
        self.wavefunction = wavefunction(i)

        #self.j_tot = math.sqrt(self.j*(self.j+1))
        
        if self.m_s < 0:
            self.spin_state = vec(0,1)
        elif self.m_s > 0:
            self.spin_state = vec(1,0)

    def data_tuple(self):
        return (self.shell, self.degeneracy, self.E_nx_ny, self.magic_number, 
                self.energy, self.nx_ny, self.n_m, self.wavefunction)

#Fill the table


s_lst = [-0.5,0.5]
#attempt to generalize fill_arrays rather than hard code
def fill_arrays(n, quantities):
    count = 0

    for i in range(n-1):
        #create list of possible quantum numbers
        #i = quantum number arrays index
        nx_ny_lst = nx_ny(i)
        n_m_lst = n_m(i+1)
        n_m_lst = n_m_lst[::-1]
            

        for j in range(i+1):
            
            for s in s_lst:
                #create instance of particle, each with unique quantum numbers
                n = n_m_lst[j].x
                m = n_m_lst[j].y
                particle = single_particle(i+1,n,m,s,count)
                particle.nx_ny = nx_ny_lst[j]
                particle.n_m = n_m_lst[j]
                #particle.m_s = k
                #particle.state = count
                particles.append(particle)
                data = (particle.shell, particle.degeneracy, particle.E_nx_ny, particle.magic_number, 
                    particle.energy, particle.nx_ny, particle.n_m, particle.wavefunction)
                count += 1

                for l in range(len(quantities)-1):
                    #store quantities in different places;
                    #for data that have type int, such as energy, E_nx_ny, etc. I put them in arrays, inserted in the 'else' block
                    #for data with other types, such as quantum number (tuple), I put them in list, inserted in the 'if' block
                    if isinstance(quantities[l],list):
                        quantities[l].append(data[l])
                    else:
                        quantities[l][i] = data[l]

#Misnomer: shells=4 actually correspond to 3 shells
shells=5
data_dict = {}

shell_arr = np.empty(shells-1)
degeneracy_arr = np.empty(shells-1)
E_nx_ny_arr = np.empty(shells-1)
magic_number_arr = np.empty(shells-1)
energy_arr = np.empty(shells-1) 
nx_ny_arr = []
n_m_arr = []
wavefunction_arr = []
particles = []
quantities = [shell_arr, degeneracy_arr, E_nx_ny_arr, magic_number_arr, energy_arr, nx_ny_arr, n_m_arr, wavefunction_arr]
names = ("Shell", "Degeneracy", "E_nx_ny", "Magic Numbers", "Energy", "(nx, ny)", "(n, m)", "Wavefunctions")


fill_arrays(shells,quantities)



'''write to files'''
'''nucleispnumbers is a file with the values: state, n, m, 1, 2*m_s'''
'''nucletwobody is a file with the values: state 1 (p), state 2 (q), spin state of p * q'''


file_spnumbers = open('nucleispnumbers.dat', 'w')
file_twobody = open('nucleitwobody.dat', 'w')
#header = "|n>,n,m,2s,2m_s"
#header = "shell,n,l,j,m,m_s"
header = "|i>,(n,m), spin state,m_s,E"
#file.write(header)
#file.write("\n")
print(header)

for p in particles:
    #q_str = "|{}>, {}, {}, {}".format(p.state,p.n_m,1,2*p.m_s)
    s_str = "|{}>, ({},{}), {}, {}, {}".format(p.state,p.n,p.m,p.spin_state, 2*p.m_s,p.shell)
    #file.write(quantum_no_str)
    #file.write("\n")
    sp_str = "{},{},{},{},{},{}".format(p.state,p.n,p.m,1,int(2*p.m_s),p.shell)
    file_spnumbers.write(sp_str)
    file_spnumbers.write("\n")
    for q in particles:
        two_body_str = "{},{},{},".format(p.state,q.state,p.spin_state.dot(q.spin_state))
        file_twobody.write(two_body_str)
        file_twobody.write("\n")
#         for r in particles:

#             for s in particles:
#                 H_0 = p.shell + q.shell + r.shell + s.shell
#                 H = H_0 + 0
#                 two_body_str = "{},{},{},{},{}".format(p.state,q.state,r.state,s.state,H)
#                 file_twobody.write(two_body_str)
#                 file_twobody.write("\n")
#     #print(q_str)
#     #print(quantum_no_str)
    print(s_str)
file_spnumbers.close()
file_twobody.close()
