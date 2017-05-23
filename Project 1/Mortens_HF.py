import numpy as np 
from decimal import Decimal
# expectation value for the one body part, Harmonic oscillator in three dimensions
def onebody(i, n, l):
    homega = 10.0
    return homega*(2*n[i] + l[i] + 1.5)

if __name__ == '__main__':

    Nparticles = 6
    """ Read quantum numbers from file """
    index = []
    n = []
    l = []
    j = []	
    mj = []
    tz = []
    spOrbitals = 0
    with open("nucleispnumbers.dat", "r") as qnumfile:
        for line in qnumfile:
            nums = line.strip('\n').split(",")
#             print(nums,type(nums),len(nums))
            if len(nums) != 0:
                index.append(int(nums[0]))
                n.append(int(nums[1]))
                l.append(int(nums[2])) #m
                j.append(int(nums[3])) #s
                mj.append(int(nums[4])) #ms
                #tz.append(int(nums[5]))
                spOrbitals += 1


        """ Read two-nucleon interaction elements (integrals) from file, brute force 4-dim array """
        nninteraction = np.zeros([spOrbitals, spOrbitals, spOrbitals, spOrbitals])
    #with open("nucleitwobody.dat","r") as infile:
    with open("coulomb.dat", "r") as infile:
        for line in infile:
            number = line.strip("\n").split(",")
            a = int(number[0])
            b = int(number[1])
            c = int(number[2])
            d = int(number[3])
            #print a, b, c, d, float(l[4])
            nninteraction[a][b][c][d] = Decimal(number[4])
    """ Set up single-particle integral """
    singleparticleH = np.zeros(spOrbitals)
    for i in range(spOrbitals):
        singleparticleH[i] = Decimal(onebody(i, n, l))

    """ Star HF-iterations, preparing variables and density matrix """

    """ Coefficients for setting up density matrix, assuming only one along the diagonals """
    C = np.eye(spOrbitals) # HF coefficients
    DensityMatrix = np.zeros([spOrbitals,spOrbitals])
    for gamma in range(spOrbitals):
        for delta in range(spOrbitals):
            sum_ = 0.0
            for i in range(Nparticles):
                sum_ += C[gamma][i]*C[delta][i]
#                 print("gamma: {} i: {} delta: {}".format(gamma,i,delta))
            DensityMatrix[gamma][delta] = Decimal(sum_)
    maxHFiter = 10
    epsilon =  1.0e-5 
    difference = 1.0
    hf_count = 0
    oldenergies = np.zeros(spOrbitals)
    newenergies = np.zeros(spOrbitals)
    while hf_count < maxHFiter and difference > epsilon:
        #print("############### Iteration %i ###############".format(hf_count))
        HFmatrix = np.zeros([spOrbitals,spOrbitals])
        for alpha in range(spOrbitals):
            for beta in range(spOrbitals):
                """  If tests for three-dimensional systems, including isospin conservation """
                if l[alpha] != l[beta] and j[alpha] != j[beta] and mj[alpha] != mj[beta] and tz[alpha] != tz[beta]: continue
                """  Setting up the Fock matrix using the density matrix and antisymmetrized NN interaction in m-scheme """
                sumFockTerm = 0.0
                for gamma in range(spOrbitals):
                    for delta in range(spOrbitals):
                        #if (mj[alpha]+mj[gamma]) != (mj[beta]+mj[delta]) and (tz[alpha]+tz[gamma]) != (tz[beta]+tz[delta]): continue
                        sumFockTerm += DensityMatrix[gamma][delta]*nninteraction[alpha][gamma][beta][delta]
                HFmatrix[alpha][beta] = Decimal(sumFockTerm)
                """  Adding the one-body term, here plain harmonic oscillator """
                if beta == alpha:   
                    HFmatrix[alpha][alpha] += singleparticleH[alpha]
                    spenergies, C = np.linalg.eigh(HFmatrix)
                """ Setting up new density matrix in m-scheme """
                DensityMatrix = np.zeros([spOrbitals,spOrbitals])
                for gamma in range(spOrbitals):
                    for delta in range(spOrbitals):
                        sum_ = 0.0
                        for i in range(Nparticles):
                            sum_ += C[gamma][i]*C[delta][i]
                        DensityMatrix[gamma][delta] = Decimal(sum_)
                newenergies = spenergies
                """ Brute force computation of difference between previous and new sp HF energies """
                sum_ =0.0
                for i in range(spOrbitals):
                    sum_ += (abs(newenergies[i]-oldenergies[i]))/spOrbitals
                difference = sum_
                oldenergies = newenergies
#                 print("Single-particle energies, ordering may have changed")
#                 for i in range(spOrbitals):
#                     print('{0:4d}  {1:.4f}'.format(i, Decimal(oldenergies[i])))
        hf_count += 1

print(spenergies)