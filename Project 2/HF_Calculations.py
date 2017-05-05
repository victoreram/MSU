from Coulomb import Coulomb
import numpy as np
import math
from decimal import Decimal


def single_H(hbar_omega, shell):
    return hbar_omega*(shell+1)

states = []
n = []
m = []
m_s = []
non_zero_combinations = []
NParticles=6
epsilon = 10**(-4)
hw=1

singleparticleH = []
with open("nucleispnumbers.dat", "r") as qnumfile:
    for line in qnumfile:
        nums = line.split(",")
        states.append(nums[0])
        n.append(int(nums[1]))
        m.append(int(nums[2]))
        m_s.append(int(nums[4]))
        singleparticleH.append(int(nums[5]))

spOrbitals = len(states)
two_interaction = np.zeros([spOrbitals,spOrbitals])

with open("nucleitwobody.dat", "r") as twobodyfile:
    for line in twobodyfile:
        nums = line.split(",")
        '''Matrix two_interaction has row, column indices for the first and second interacting electron respectively'''
        '''Value at those indices is the inner product of the two spin states'''
        two_interaction[int(nums[0])][int(nums[1])] = int(nums[2])

hbar_omega = 1
CMatrix = np.eye(spOrbitals) # HF coefficients
DensityMatrix = np.zeros([spOrbitals,spOrbitals])
HFMatrix = np.zeros([spOrbitals,spOrbitals])
Coulomb.write_coulomb_file(hw,nuclei_file="nucleispnumbers.dat")
CoulombMatrix = Coulomb.read_coulomb_file(spOrbitals,coulomb_file="coulomb.dat")
for gamma in range(spOrbitals):
    for delta in range(spOrbitals):
        sum_ = 0.0
        for i in range(NParticles):
            sum_ += CMatrix[gamma][i]*CMatrix[delta][i]
#                 print("gamma: {} i: {} delta: {}".format(gamma,i,delta))
        DensityMatrix[gamma][delta] = Decimal(sum_)

print(DensityMatrix)

hf_count = 0
maxHFiter = 10
difference = 1
epsilon = 10**-3
# for i in range(spOrbitals):
#     '''Fix this: 1,1,2,2,2,2,etc.'''
#     #singleparticleH[i] = Decimal(single_H(hbar_omega,i))
#     singleparticleH[i] = 1
    
with open("hf_energies.txt", "w") as hffile:
    oldenergies = np.zeros(spOrbitals)
    newenergies = np.zeros(spOrbitals)
    while hf_count < maxHFiter and difference > epsilon:
        HFMatrix = np.zeros([spOrbitals,spOrbitals])
        for alpha in range(spOrbitals):
            for beta in range(spOrbitals):
                
                M_s_ab = m_s[alpha] + m_s[beta]
                M_ab = m[alpha] + m[beta]
                '''Add initial term for E_a_b'''
                if beta == alpha:   
                    HFMatrix[alpha][beta] += singleparticleH[alpha]
                sum_1 = 0.0
                for gamma in range(spOrbitals):
                    for delta in range(spOrbitals):
                        M_s_cd = m_s[gamma] + m_s[delta]
                        M_cd = m[gamma] + m[delta]
                        C_sum = 0.0
                        direct_exchange_term = 0.0

                        '''Test for spin and M conservation'''
                        if M_s_ab == M_s_cd and M_ab == M_cd:
                            direct = two_interaction[alpha][gamma]*two_interaction[beta][delta]*CoulombMatrix[alpha][beta][gamma][delta]
                            exchange = two_interaction[alpha][delta]*two_interaction[beta][gamma]*CoulombMatrix[alpha][beta][delta][gamma]

                            '''Direct *- Coulomb(alpha,beta,gamma,delta)'''
                            direct_exchange_term = (direct - exchange)
                            
                            #print(gamma,delta,direct_exchange_term)
                            #print(gamma,delta,DensityMatrix[gamma][delta])
                            sum_1 += DensityMatrix[gamma][delta]*direct_exchange_term       
                HFMatrix[alpha][beta] += sum_1                 

                #print(sum_1)

        spenergies, CMatrix = np.linalg.eigh(HFMatrix)
        #print(spenergies)
        for gamma in range(spOrbitals):
            for delta in range(spOrbitals):
                sum_ = 0.0
                for i in range(NParticles):
                    sum_ += CMatrix[gamma][i]*CMatrix[delta][i]
        #                 print("gamma: {} i: {} delta: {}".format(gamma,i,delta))
                DensityMatrix[gamma][delta] = Decimal(sum_)
#                 '''Summing C terms'''
#                 for j in range(NParticles):
#                     C_sum += CMatrix[j][gamma]*CMatrix[j][delta]

#                 '''Update Density and HF Matrix'''
#                 DensityMatrix[gamma][delta] = Decimal(C_sum)
    #spenergies, CMatrix = np.linalg.eigh(HFMatrix)
        newenergies = spenergies
        """ Brute force computation of difference between previous and new sp HF energies """
        sum_ =0.0
        for i in range(spOrbitals):
            sum_ += (abs(newenergies[i]-oldenergies[i]))/spOrbitals
        difference = sum_
        oldenergies = newenergies

        #print("difference ", difference)

        hf_count += 1
        #print(HFMatrix)
    #hffile.write(newenergies)
#print("Final Density Matrix \n {}".format(DensityMatrix))
        print("SP energies \n {}".format(spenergies))
    print(DensityMatrix)
