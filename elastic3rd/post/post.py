#!python
#

import numpy as np
import esfit

def get_cij(coef_fit, crystaltype):
    (coef2, coef3) = escoef(crystaltype)
    C3 = np.linalg.solve(coef3, 6.0*coef_fit[:, 1])
    #C2 = np.linalg.solve(coef2[[0, 1, 6], [0, 1, 6]], 2.0*coef_fit[[0, 1, 6], 0])
    C2 = np.linalg.lstsq(coef2, 2.0*coef_fit[:, 0], rcond = None)[0]
    return (C2, C3)

def get_coef(s, e, V0):
    #s: strain, e:energy
    eVpmol2GPa = 160.21719175
    (m, n) = e.shape
    coef_fit = np.zeros((n, 2))
    for i in range(0, n):
        e0 = e[int((m-1)/2)][i]
        for j in range(0, m):
            e[j][i] = (e[j][i] - e0)/V0*eVpmol2GPa
        (coefi, pcovi) = esfit.esfit(s, e[:, i])
        coef_fit[i, :] = coefi
    return coef_fit

def read_e(EEnergy = "EEnergy.txt"):
    e = np.loadtxt(EEnergy)
    return e

def escoef(crystaltryp):
    if crystaltryp == "cubic":
        #REF: Phys. Rev. B, 75(2007), 094105
        #coef2 order: C11, C12, C44
        coef2 = np.array([[1, 0, 0], [2, 2, 0], [3, 6, 0], 
            [1, 0, 4], [1, 0, 4], [0, 0, 12]])
        #coef3 order: C111, C112, C123, C144, C166, C456
        coef3 = np.array([[1, 0, 0, 0, 0, 0], [2, 6, 0, 0, 0, 0], 
            [3, 18, 6, 0, 0, 0], [1, 0, 0, 12, 0, 0], 
            [1, 0, 0, 0, 12, 0], [0, 0, 0, 0, 0, 48]])
    elif crystaltryp == "trigonal":
        #ref: Phys. Rev. B, 75(2007), 094105
        #coef2 order: C11, C12, C13, C14, C33, C44
        coef2 = np.array([[1, 0, 0, 0, 0, 0], [2, 2, 0, 0, 0, 0], 
            [2, 2, 4, 0, 1, 0], [1, 0, 0, 4, 0, 4], [1, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 1, 0], [1, 0, 2, 0, 1, 0], [0, 0, 0, 0, 0, 4],
            [2, -2, 0, 0, 1, 0], [1, 0, 0, -4, 0, 4], [1, 0, 0, 0, 0, 4],
            [0, 0, 0, 0, 1, 4], [1, 0, 2, 4, 1, 4], [2, 2, 0, 0, 0, 4]])
        #coef3 order: C111, C112, C113, C114, C123, C124, C133
        #             C134, C144, C155, C222, C333, C344, C444
        coef3 = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [4, 6, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0], 
            [4, 6, 6, 0, 6, 0, 6, 0, 0, 0, -2, 1, 0, 0], 
            [1, 0, 0, 6, 0, 0, 0, 0, 6, 0, 0, 0, 0, 8], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], 
            [0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8], 
            [0, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 1, 0, 0], 
            [0, 0, 0, -6, 0, -12, 0, 0, 0, 12, 1, 0, 0, 8], 
            [0, 0, 0, 0, 0, 0, 0, 0, 12, 0, 1, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 12, 0], 
            [1, 0, 3, 6, 0, 0, 3, 12, 12, 0, 0, 1, 12, 8], 
            [4, 6, 0, 0, 0, 0, 0, 0, 12, 12, -2, 0, 0, 0]])
    return (coef2, coef3)

'''
crystaltype = "cubic"

e = read_e()
print e
s = np.arange(-0.05, 0.075, 0.025)
print s
V0 = 163.125598397
coef_fit = get_coef(s, e, V0)
print coef_fit

(C2, C3) = get_cij(coef_fit, crystaltype)
print C2
print C3
'''