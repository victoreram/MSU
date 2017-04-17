# #include "Coulomb_Functions.h"
import math
import numpy as np


def logfac(n):
    fac = 0.0
    for i in range(2,n+1):
        fac += np.log(i)
    return fac

def logratio1(n1,n2,n3,n4):
    return -logfac(n1) - logfac(n2) - logfac(n3) - logfac(n4)

def logratio2(G):
    return -0.5*(G+1)*math.log(2)

def product1(n1,n2,n3,n4,m1,m2,m3,m4):
    prod = logfac(n1) + logfac(n2) + logfac(n3) + logfac(n4)
    arg1 = n1 + abs(m1)
    arg2 = n2 + abs(m2)
    arg3 = n3 + abs(m3)
    arg4 = n4 + abs(m4)
    prod -= (logfac(arg1)+ logfac(arg2) + logfac(arg3)+ logfac(arg4))
    prod *= 0.5
    return math.exp(prod)

def logproduct2(n1,n2,n3,n4,m1,m2,m3,m4,j1,j2,j3,j4):
    arg1 = n1 + abs(m1)
    arg2 = n2 + abs(m2)
    arg3 = n3 + abs(m3)
    arg4 = n4 + abs(m4)
    narg1 = n1 - j1
    narg2 = n2 - j2
    narg3 = n3 - j3
    narg4 = n4 - j4
    jarg1 = j1 + abs(m1)
    jarg2 = n2 + abs(m2)
    jarg3 = n3 + abs(m3)
    jarg4 = n4 + abs(m4)
    prod = logfac(arg1) + logfac(arg2) + logfac(arg3) + logfac(arg4)
    prod -= (logfac(narg1) + logfac(narg2) + logfac(narg3) + logfac(narg4));
    prod -= (logfac(jarg1) + logfac(jarg2) + logfac(jarg3) + logfac(jarg4));
    return prod

def logproduct3(l1, l2, l3, l4, g1, g2, g3, g4):
    garg1 = g1 - l1
    garg2 = g2 - l2
    garg3 = g3 - l3
    garg4 = g4 - l4
    prod = logfac(g1) + logfac(g2) + logfac(g3) + logfac(g4)
    prod -= (logfac(l1) + logfac(l2) + logfac(l3) + logfac(l4))
    prod -= (logfac(garg1) + logfac(garg2) + logfac(garg3) + logfac(garg4))
    return prod

def Coulomb_HO(hw, ni, mi, nj, mj, nk, mk, nl, ml):
    direct = 0.0
    exchange = 0.0
    if((mi + mj) != (mk + ml)):
        return 0.0
    for j1 in range(0,ni+1):
        for j2 in range(0,nj+1):
            for j3 in range(0, nl+1):
                for j4 in range(0, nk+1):
                    g1 = int(j1 + j4 + 0.5*(abs(mi) + mi) + 0.5*(abs(mk) - mk))
                    g2 = int(j2 + j3 + 0.5*(abs(mj) + mj) + 0.5*(abs(ml) - ml))
                    g3 = int(j3 + j2 + 0.5*(abs(ml) + ml) + 0.5*(abs(mj) - mj))
                    g4 = int(j4 + j1 + 0.5*(abs(mk) + mk) + 0.5*(abs(mi) - mi))
                    G = g1 + g2 + g3 + g4
                    LogRatio1 = logratio1(j1, j2, j3, j4)
                    LogProd2 = logproduct2(ni, mi, nj, mj, nl, ml, nk, mk, j1, j2, j3, j4)
                    LogRatio2 = logratio2(G)
                    temp = 0.0
                    for l1 in range(0,g1+1):
                        for l2 in range(0,g2+1):
                            for l3 in range(0, g3+1):
                                for l4 in range(0, g4+1):
                                    if((l1 + l2) != (l3 + l4)):
                                        continue
                                    L = l1 + l2 + l3 + l4
                                    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1) * math.exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + math.lgamma(1.0 + 0.5*L) + math.lgamma(0.5*(G - L + 1.0)))
                    direct += (-2*((j1 + j2 + j3 + j4)%2) + 1) * math.exp(LogRatio1 + LogProd2 + LogRatio2) * temp

    direct *= product1(ni, mi, nj, mj, nl, ml, nk, mk)
    return math.sqrt(hw)*(direct - exchange)

def write_coulomb_file(hw=1,nuclei_file="nucleispnumbers.dat"):

    states = []
    n = []
    m = []
    m_s = []
    with open(nuclei_file, "r") as qnumfile:
        for line in qnumfile:
            nums = line.split(",")
            states.append(nums[0])
            n.append(int(nums[1]))
            m.append(int(nums[2]))
            m_s.append(int(nums[4]))
            
    spOrbitals = len(states)
    coulomb_file = open('coulomb.dat', 'w')
    V = np.zeros((spOrbitals,spOrbitals,spOrbitals,spOrbitals))
    for i in range(spOrbitals):
        ni = n[i]
        mi = m[i]
        for j in range(spOrbitals):
            nj = n[j]
            mj = m[j]
            for k in range(spOrbitals):
                nk = n[k]
                mk = m[k]
                for l in range(spOrbitals):
                    nl = n[l]
                    ml = m[l]
                    V_ijkl = Coulomb_HO(hw,ni,mi,nj,mj,nk,mk,nl,ml)
                    V[i][j][k][l] = V_ijkl
                    coulomb_str = "{},{},{},{},{},".format(str(i),str(j),str(k),str(l),str(V_ijkl))
                    coulomb_file.write(coulomb_str)
                    coulomb_file.write("\n")
                    

    return V
                    
def read_coulomb_file(spOrbitals,coulomb_file="coulomb.dat"):
    '''Reads coulomb file, returns a matrix of potentials'''
    read_file = open(coulomb_file,"r")
    VMatrix = np.zeros((spOrbitals,spOrbitals,spOrbitals,spOrbitals),dtype=np.float)
    for line in read_file:
        nums = line.split(",")
        i = int(nums[0])
        j = int(nums[1])
        k = int(nums[2])
        l = int(nums[3])
        V = float(nums[4])
        VMatrix[i][j][k][l] = V
        
    return VMatrix
                    
