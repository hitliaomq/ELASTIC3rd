#!python
#filename: symdata.py
#function: the independent elastic constant for different symmeties
#author: Mingqing Liao
#e-mail: liaomq1900127@163.com
#FGMS group @ Harbin Institute of Technology(HIT)

import numpy as np

def coef_crystal(CrystalType = "c1", Ord = 3):
    CrystalType = CrystalType.lower()
    if Ord == 2:
        A = [0.5,-0.5]
        if CrystalType == "triclinic" or CrystalType == "n":
            CoefCoef = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        elif CrystalType == "monoclinic" or CrystalType == "m":
            CoefCoef = [1,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,1,0,1]
        elif CrystalType == "orthorhombic" or CrystalType == "o":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1]
        elif CrystalType == "tetragonal2" or CrystalType == "t2":
            CoefCoef = [1,1,1,0,0,1,1,1,0,0,-1,1,0,0,0,1,0,0,1,0,1]
        elif CrystalType == "tetragonal1" or CrystalType == "t1":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1]        
        elif CrystalType == "rhombohedral2" or CrystalType == "r2":
            CoefCoef = [1,1,1,1,1,0,1,1,-1,-1,0,1,0,0,0,1,0,-1,1,1,A]
        elif CrystalType == "rhombohedral1" or CrystalType == "r1":
            CoefCoef = [1,1,1,1,0,0,1,1,-1,0,0,1,0,0,0,1,0,0,1,1,A]
        elif CrystalType == "hexagonal2" or CrystalType == "h2":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,A]
        elif CrystalType == "hexagonal1" or CrystalType == "h1":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,A]        
        elif CrystalType == "cubic2" or CrystalType == "c2":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1]
        elif CrystalType == "cubic1" or CrystalType == "c1":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1]
        elif CrystalType == "isotropic" or CrystalType == "i":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,A,0,0,A,0,A]
    elif Ord == 3:
        A = [1,1,-1];B = [-0.5,-1.5];C = [0.5,1.5];D = [-0.5,-0.25,0.75];E = [-1,-2];F = [-1,-2]
        G = [-0.5,0.5];H = [0.5,-0.5];I = [0.5,-0.25,-0.25];J = [0.5,-0.5];K = [-0.5,0.5]
        L = [0.5,-0.5];M = [0.25,-0.25];N = [0.125,-0.375,0.25]
        if CrystalType == "triclinic" or CrystalType == "n":
            CoefCoef = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        elif CrystalType == "monoclinic" or CrystalType == "m":
            CoefCoef = [1,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,
                1,0,1,1,0,1,1,0,1,0,1,0,1,1,0,1,0,1,0,0,1,0,1,0,1,0]
        elif CrystalType == "orthorhombic" or CrystalType == "o":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,0,
                1,0,0,1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0]
        elif CrystalType == "tetragonal2" or CrystalType == "t2":
            CoefCoef = [1,1,1,0,0,1,1,1,0,0,0,1,0,0,1,1,1,0,1,0,1,1,1,0,0,-1,1,0,0,-1,
                1,-1,0,1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,-1,0,0]
        elif CrystalType == "tetragonal1" or CrystalType == "t1":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,0,1,
                0,0,1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0]        
        elif CrystalType == "rhombohedral2" or CrystalType == "r2":
            CoefCoef = [1,1,1,1,1,1,A,1,1,1,-1,1,1,1,0,1,1,B,1,C,D,1,1,E,F,1,1,-1,-1,0,
                1,-1,G,1,H,I,1,0,0,0,1,0,-1,1,1,J,1,1,1,-1,K,1,-1,-1,1,-1]
        elif CrystalType == "rhombohedral1" or CrystalType == "r1":
            CoefCoef = [1,1,1,1,0,0,A,1,1,0,0,1,1,0,0,1,0,0,1,C,D,1,1,E,0,0,1,-1,0,0,1,
                0,0,1,H,I,1,0,0,0,1,0,0,1,1,J,1,0,0,-1,K,1,0,0,0,0]
        elif CrystalType == "hexagonal2" or CrystalType == "h2":
            CoefCoef = [1,1,1,0,0,1,A,1,0,0,-1,1,0,0,0,1,1,0,1,0,D,1,1,0,0,1,1,0,0,0,1,
                -1,0,1,0,I,1,0,0,0,1,0,0,1,0,J,0,0,1,0,K,0,0,-1,0,-1]
        elif CrystalType == "hexagonal1" or CrystalType == "h1":
            CoefCoef = [1,1,1,0,0,0,A,1,0,0,0,1,0,0,0,1,0,0,1,0,D,1,1,0,0,0,1,0,0,0,1,
                0,0,1,0,I,1,0,0,0,1,0,0,1,0,J,0,0,0,0,K,0,0,0,0,0]        
        elif CrystalType == "cubic2" or CrystalType == "c2":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,0,1,
                0,0,1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0]
        elif CrystalType == "cubic1" or CrystalType == "c1":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,0,1,
                0,0,1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0]
        elif CrystalType == "isotropic" or CrystalType == "i":
            CoefCoef = [1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,L,0,0,M,0,M,1,1,0,0,0,1,0,0,0,M,
                0,0,L,0,M,1,0,0,0,M,0,0,M,0,L,0,0,0,0,N,0,0,0,0,0]
    return CoefCoef

def coef_ind(CrystalType = "c1", Ord = 3):
    #ATTENTION: here the first element is 1, not 0
    #for 2nd elastic constant
    CrystalType = CrystalType.lower()
    if Ord == 2:
        A = [1, 2]
        if CrystalType == "triclinic" or CrystalType == "n":
            CoefInd = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
        elif CrystalType == "monoclinic" or CrystalType == "m":
            CoefInd = [1,2,3,0,5,0,7,8,0,10,0,12,0,14,0,16,0,18,19,0,21]
        elif CrystalType == "orthorhombic" or CrystalType == "o":
            CoefInd = [1,2,3,0,0,0,7,8,0,0,0,12,0,0,0,16,0,0,19,0,21]
        elif CrystalType == "tetragonal2" or CrystalType == "t2":
            CoefInd = [1,2,3,0,0,6,1,3,0,0,6,12,0,0,0,16,0,0,16,0,21]
        elif CrystalType == "tetragonal1" or CrystalType == "t1":
            CoefInd = [1,2,3,0,0,0,1,3,0,0,0,12,0,0,0,16,0,0,16,0,21]        
        elif CrystalType == "rhombohedral2" or CrystalType == "r2":
            CoefInd = [1,2,3,4,5,0,1,3,4,5,0,12,0,0,0,16,0,5,16,4,A]
        elif CrystalType == "rhombohedral1" or CrystalType == "r1":
            CoefInd = [1,2,3,4,0,0,1,3,4,0,0,12,0,0,0,16,0,0,16,4,A]
        elif CrystalType == "hexagonal2" or CrystalType == "h2":
            CoefInd = [1,2,3,0,0,0,1,3,0,0,0,12,0,0,0,16,0,0,16,0,A]
        elif CrystalType == "hexagonal1" or CrystalType == "h1":
            CoefInd = [1,2,3,0,0,0,1,3,0,0,0,12,0,0,0,16,0,0,16,0,A]
        elif CrystalType == "cubic2" or CrystalType == "c2":
            CoefInd = [1,2,2,0,0,0,1,2,0,0,0,1,0,0,0,16,0,0,16,0,16]
        elif CrystalType == "cubic1" or CrystalType == "c1":
            CoefInd = [1,2,2,0,0,0,1,2,0,0,0,1,0,0,0,16,0,0,16,0,16]
        elif CrystalType == "isotropic" or CrystalType == "i":
            CoefInd = [1,2,2,0,0,0,1,2,0,0,0,1,0,0,0,A,0,0,A,0,A]
    elif Ord == 3:
        A = [1, 2, 22];B = [5, 1];C = [4, 9];D = [1, 2, 22];E = [4, 9];F = [5, 10];G = [5, 10]
        H = [4, 9];I = [1, 2, 22];J = [3, 8];K = [16, 19];L = [2, 8,];M = [1, 2];N = [1, 2, 8]
        if CrystalType == "triclinic" or CrystalType == "n":
            CoefInd = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
                29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56]
        elif CrystalType == "monoclinic" or CrystalType == "m":
            CoefInd = [1,2,3,0,5,0,7,8,0,10,0,12,0,14,0,16,0,18,19,0,21,22,23,0,25,0,27,
                0,29,0,31,0,33,34,0,36,37,0,39,0,41,0,43,44,0,46,0,48,0,0,51,0,53,0,55,0]
        elif CrystalType == "orthorhombic" or CrystalType == "o":
            CoefInd = [1,2,3,0,0,0,7,8,0,0,0,12,0,0,0,16,0,0,19,0,21,22,23,0,0,0,27,0,0,
                0,31,0,0,34,0,36,37,0,0,0,41,0,0,44,0,46,0,0,0,0,51,0,0,0,0,0]
        elif CrystalType == "tetragonal2" or CrystalType == "t2":
            CoefInd = [1,2,3,0,0,6,2,8,0,0,0,12,0,0,15,16,17,0,19,0,21,1,3,0,0,6,12,0,0,
                15,19,17,0,16,0,21,37,0,0,0,41,0,0,41,0,46,0,0,49,0,51,0,0,49,0,0]
        elif CrystalType == "tetragonal1" or CrystalType == "t1":
            CoefInd = [1,2,3,0,0,0,2,8,0,0,0,12,0,0,0,16,0,0,19,0,21,1,3,0,0,0,12,0,0,0,
                19,0,0,16,0,21,37,0,0,0,41,0,0,41,0,46,0,0,0,0,51,0,0,0,0,0]        
        elif CrystalType == "rhombohedral2" or CrystalType == "r2":
            CoefInd = [1,2,3,4,5,6,A,8,9,10,6,12,13,14,0,16,17,B,19,C,D,22,3,E,F,6,12,13,
                14,0,19,17,G,16,H,I,37,0,0,0,41,0,14,41,13,J,47,48,17,47,K,9,48,17,10,6]
        elif CrystalType == "rhombohedral1" or CrystalType == "r1":
            CoefInd = [1,2,3,4,0,0,A,8,9,0,0,12,13,0,0,16,0,0,19,C,D,22,3,E,0,0,12,13,0,0,
                19,0,0,16,H,I,37,0,0,0,41,0,0,41,13,J,47,0,0,47,K,9,0,0,0,0]
        elif CrystalType == "hexagonal2" or CrystalType == "h2":
            CoefInd = [1,2,3,0,0,6,A,8,0,0,6,12,0,0,0,16,17,0,19,0,D,22,3,0,0,6,12,0,0,0,
                19,17,0,16,0,I,37,0,0,0,41,0,0,41,0,J,0,0,17,0,K,0,0,17,0,6]
        elif CrystalType == "hexagonal1" or CrystalType == "h1":
            CoefInd = [1,2,3,0,0,0,A,8,0,0,0,12,0,0,0,16,0,0,19,0,D,22,3,0,0,0,12,0,0,0,
                19,0,0,16,0,I,37,0,0,0,41,0,0,41,0,J,0,0,0,0,K,0,0,0,0,0]
        elif CrystalType == "cubic2" or CrystalType == "c2":
            CoefInd = [1,2,3,0,0,0,3,8,0,0,0,2,0,0,0,16,0,0,19,0,21,1,2,0,0,0,3,0,0,0,21,
                0,0,16,0,19,1,0,0,0,19,0,0,21,0,16,0,0,0,0,51,0,0,0,0,0]
        elif CrystalType == "cubic1" or CrystalType == "c1":
            CoefInd = [1,2,2,0,0,0,2,8,0,0,0,2,0,0,0,16,0,0,19,0,19,1,2,0,0,0,2,0,0,0,19,
                0,0,16,0,19,1,0,0,0,19,0,0,19,0,16,0,0,0,0,51,0,0,0,0,0]
        elif CrystalType == "isotropic" or CrystalType == "i":
            CoefInd = [1,2,2,0,0,0,2,8,0,0,0,2,0,0,0,L,0,0,M,0,M,1,2,0,0,0,2,0,0,0,M,0,0,
                L,0,M,1,0,0,0,M,0,0,M,0,L,0,0,0,0,N,0,0,0,0,0]
    return CoefInd

def get_unique(CoefInd):
    #ATTENTION: here the first element is 1, not 0
    n = len(CoefInd)
    CoefUniq = []
    for i in range(0, n):
        coefi = CoefInd[i]
        if type(coefi) is int:
            if coefi == 0:
                pass
            else:
                #coefi = coefi - 1
                if not coefi in CoefUniq:
                    CoefUniq.append(coefi)
    #print len(CoefUniq)
    return CoefUniq

def group2crytyp(group):
    if group == "1" or group == "-1":
        CrystalType = "triclinic"
    elif group == "2" or group == "m" or group == "2/m":
        CrystalType = "monoclinic"
    elif group == "222" or group == "mm2" or group == "222/mmm":
        CrystalType = "orthorhombic"
    elif group == "4" or group == "-4" or group == "4/m":
        CrystalType = "tetragonal2"
    elif group == "422" or group == "4mm" or group == "-42m" or group == "422/mmm":
        CrystalType = "tetragonal1"
    elif group == "23" or group == "2/m-3":
        CrystalType = "cubic2"
    elif group == "432" or group == "-43m" or group == "4/m-32/m":
        CrystalType = "cubic1"
    elif group == "3" or group == "-3":
        CrystalType = "rhombohedral2"
    elif group == "32" or group == "3m" or group == "-32/m":
        CrystalType = "rhombohedral1"
    elif group == "6" or group == "-6" or group == "6/m":
        CrystalType = "hexagonal2"
    elif group == "622" or group == "6mm" or group == "-6m2" or group == "622/mmm":
        CrystalType = "hexagonal1"
    return CrystalType

'''
CrystalType = "h1"
CoefInd = coef_ind2(CrystalType)
CoefUniq = get_unique(CoefInd)
print CoefUniq
print len(CoefUniq)
#UniqeCoef = np.unique(IndepCoef)
print CoefInd
#print UniqeCoef
print len(CoefInd)

a = [1, [1, 2, 3], 2, 3, 4]

print type(a[1]) is list
'''

'''
CrystalType = "r1"
Ord = 4
coefcoef4 = coef_crystal(CrystalType, Ord)
print coefcoef4
print len(coefcoef4)
coefind4 = coef_ind(CrystalType, Ord)
print coefind4
print len(coefind4)
'''
