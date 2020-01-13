#!python
#filename: symmetry.py
#function: generate the strain patterns 
#          for calculating third-order elastic constant using strain-energy method
#author: Mingqing Liao
#e-mail: liaomq1900127@163.com
#FGMS group @ Harbin Institute of Technology(HIT)

from itertools import combinations, product
import elastic3rd.crystal.deform as crydef
import numpy as np
import math
import elastic3rd.symmetry.symdata as symdata

def gen_strain(strain_sca):
    #function: generate the different strain patterns, containing 0, strain_sca and strain_sca/2
    #          if you want to change the strain patterns, modify the strainprod in line 19
    strain_sca = float(strain_sca)
    Strain = np.array([[0, 0, 0, 0, 0, 0]])
    for i in range(1, 7):
        straini = list(combinations([0, 1, 2, 3, 4, 5], i))
        for j in range(0, len(straini)):
            straintmp = np.array([[0, 0, 0, 0, 0, 0]], dtype = float)
            strainindex = np.array(straini[j], dtype = int)
            #strainprod = list(product([strain_sca], repeat=i))
            strainprod = list(product([strain_sca, strain_sca/2.], repeat=i))
            for k in range(0, len(strainprod)):
                strainprodk = np.array(strainprod[k])
                #print strainprodk
                straintmp[0, strainindex] = strainprodk         
                straintmp[0, 3:] = 2 * straintmp[0, 3:]
                Strain = np.append(Strain, straintmp, axis = 0)
    #Strain = np.append(Strain, np.array([[1, 0.5, 0, 0, 0, 0]]), axis = 0)
    Strain = np.delete(Strain, 0, 0)
    return Strain

def Num2IJK(Num, N, Ord):
    IJK = np.zeros(int(Ord), dtype=int)
    for i in range(int(Ord), 0, -1):
        k = i - 1
        DivNum = N ** k
        IJK[int(Ord) - 1 - k] = math.ceil(float(Num) / float(DivNum))
        if IJK[int(Ord) - 1 -k] > N:
            raise IndexError('ERROR')
        #print Ord - 1 - k
        Num = Num % DivNum
        if Num == 0:
            #IJK[Ord - 1 - k] = IJK[Ord - 1 - k] - 1
            for j in range(int(Ord) - k, int(Ord)):
                #print j
                IJK[j] = N
            return IJK
    return IJK

def gen_strain_coef(StrainOrd, StrainOrdCoef, Ord = 3):
    #function: get the non-zero second- or third-(specified by the Ord parameter) order 
    #          elastic constant(CoefUniq) and corresponding coefficient(CoefUniqCoef)
    #          the input parameter can get by gen_strainord
    #          Ord the order, 2..9, no more than 9
    n_strain = len(StrainOrd)
    N = int(n_strain ** Ord)
    #flag = 0
    Coef = np.ones(N)
    Cijk = np.zeros(N)
    for i in range(0, N):
        IJK = Num2IJK(i + 1, n_strain, Ord)
        #print IJK
        coeftmp = np.zeros(int(Ord))
        for i_ord in range(0, int(Ord)):
            #print IJK[i_ord]
            coeftmp[i_ord] = StrainOrd[IJK[i_ord] - 1]
        coeftmp = np.sort(coeftmp)
        #print coeftmp
        for j in range(0, int(Ord)):
            Cijk[i] = Cijk[i] + coeftmp[j] * (10 ** (Ord - 1 - j))
            Coef[i] = Coef[i] * StrainOrdCoef[IJK[j] - 1]
    CoefUniq, CoefIndex = np.unique(Cijk, return_index=True)
    #print len(CoefUniq)
    Coef = Coef[CoefIndex]
    CoefUniqCoef = np.zeros(len(CoefUniq))
    for i in range(0, len(CoefUniq)):
        CoefUniqCoef[i] = sum(Cijk == CoefUniq[i])
    CoefUniqCoef = CoefUniqCoef * Coef
    return (CoefUniq, CoefUniqCoef)

def gen_strainord(strain):
    #function: get the non-zero elements in the strain tensor(StrainOrd)
    #          and get the corresponding value(StrainOrdCoef)
    StrainOrd = []
    StrainOrdCoef = []
    if crydef.is_strain_matrix(strain):
        for i in range(0, 3):
            for j in range(0, 3):
                if not (strain[i, j] == 0):
                    StrainOrd.append(voigtmap(i+1, j+1))
                    StrainOrdCoef.append(float(strain[i][j]))
    else:
        raise IOError("The input strain matrix is not symmetry.")
    return (StrainOrd, StrainOrdCoef)


def voigtmap(i, j):
    #function: convert i,j pair to single value using Voigt's notation
    if i * j == 1:
        voigt = 1
    elif i * j == 4:
        voigt = 2
    elif i * j == 9:
        voigt = 3
    elif i * j == 6:
        voigt = 4
    elif i * j == 3:
        voigt = 5
    elif i * j == 2:
        voigt = 6
    else:
        voigt = 0
        IOError("Both i and j should range from 1 to 3")
        print("ERROR: Both i and j shold range from 1 to 3.")
    return voigt

def gen_cijk_coef(CrystalType, Strain, Ord):
    #function: generate the coefficients of independent elastic constant according to 
    #          the crystal type(CrystalType) and the strain(Strain) and order(Ord)
    (StrainOrd, StrainOrdCoef) = gen_strainord(Strain)
    CijkInd = symdata.coef_ind(CrystalType, Ord)
    CoefCoef = symdata.coef_crystal(CrystalType, Ord)                
    (Cijk, CijkCoef) = gen_strain_coef(StrainOrd, StrainOrdCoef, Ord)
    CijkUniq = np.array(symdata.get_unique(CijkInd))
    CoefStrain = np.zeros(len(CijkUniq))
    CijkNum = cijk2num(Cijk, Ord)
    for i in range(0, len(CijkNum)):
        COrderi = CijkInd[CijkNum[i]]
        #print COrderi
        #print COrderi is list
        if type(COrderi) is list:
            for j in range(0, len(COrderi)):
                #print CijkNum
                #print CoefCoef[CijkNum[i]][j]
                #print CijkCoef[i]
                CoefOrderi = float(CoefCoef[CijkNum[i]][j]) * float(CijkCoef[i])
                Index_CoefStrain = np.argwhere(CijkUniq == COrderi[j])
                Index_CoefStrain = Index_CoefStrain[0][0]
                CoefStrain[Index_CoefStrain] = CoefStrain[Index_CoefStrain] + CoefOrderi
        else:
            if COrderi != 0:
                CoefOrderi = float(CoefCoef[CijkNum[i]]) * float(CijkCoef[i])
                Index_CoefStrain = np.argwhere(CijkUniq == COrderi)
                #print Index_CoefStrain
                Index_CoefStrain = Index_CoefStrain[0][0]
                CoefStrain[Index_CoefStrain] = CoefStrain[Index_CoefStrain] + CoefOrderi
    return CoefStrain

def gen_strain_mode(CrystalType, Ord):
    #function: generate the strain mode for given crystal type(CrystalType) and order(Ord)
    StrainAll = gen_strain(1)
    CijkInd = symdata.coef_ind(CrystalType, Ord)
    CoefCoef = symdata.coef_crystal(CrystalType, Ord)
    CijkUniq = np.array(symdata.get_unique(CijkInd))
    n_uniq = len(CijkUniq)
    Cijk = num2cijk(CijkUniq, Ord)
    StrainMode = np.zeros((n_uniq, 3, 3), dtype = float)
    #StrainModeCoef = np.zeros((n_uniq, n_uniq), dtype = float)

    StrainModeCoef = CoefStr(Ord)
    exec('StrainModeCoef.coef' + str(int(Ord)) + ' = np.zeros((n_uniq, n_uniq), dtype = float)')
    for i in range(2, int(Ord)):
        CijUniq = np.array(symdata.get_unique(symdata.coef_ind(CrystalType, i)))
        n_uniqi = len(CijUniq)
        exec('StrainModeCoef.coef' + str(i) + ' = np.zeros((n_uniq, n_uniqi), dtype = float)')
    count = 0
    for i in range(0, len(StrainAll)):
        Straini = crydef.vec2matrix(StrainAll[i])
        StrainModei = gen_cijk_coef(CrystalType, Straini, Ord)
        if i == 0:
            StrainMode[count, :, :] = Straini
            exec('StrainModeCoef.coef' + str(int(Ord)) + '[count, :] = StrainModei')
            for j in range(2, int(Ord)):
                exec('StrainModeCoef.coef' + str(j) + '[count, :] = gen_cijk_coef(CrystalType, Straini, j)')               
            count = count + 1
        else:
            StrainModetmp = np.zeros((count + 1, n_uniq), dtype = float)
            #print StrainMode[0, :]
            StrainModetmp[0:count, :] = eval('StrainModeCoef.coef' + str(int(Ord)) + '[0:count, :]')
            StrainModetmp[count, :] = StrainModei
            SMRank = np.linalg.matrix_rank(StrainModetmp)
            if SMRank == count + 1:
                #print SMRank
                exec('StrainModeCoef.coef' + str(int(Ord)) + '[count, :] = StrainModei')
                StrainMode[count, :, :] = Straini
                for j in range(2, int(Ord)):
                    exec('StrainModeCoef.coef' + str(j) + '[count, :] = gen_cijk_coef(CrystalType, Straini, j)') 
                #print StrainModei
                count = count + 1
                if count == n_uniq:
                    break
    #print StrainModeCoef
    #print StrainMode
    return (StrainModeCoef, StrainMode)

def gen_cijall(Ord):
    N = int(6 ** Ord)
    CijkAll = np.zeros(N, dtype=int)
    for i in range(0, N):
        IJK = Num2IJK(i + 1, 6, int(Ord))
        IJK = np.sort(IJK)
        for j in range(0, int(Ord)):
            CijkAll[i] = CijkAll[i] + IJK[j] * (10 ** (Ord - 1 - j))
    CijkAll = np.unique(CijkAll)
    return CijkAll

def cijk2num(Cijk, Ord):
    #function: convert the cijk to number
    #ATTENTION: here the first element is 0, not 1
    CijkAll = gen_cijall(Ord)
    CijkNum = []
    for i in Cijk:
        indexi = np.argwhere(CijkAll == i)
        CijkNum.append(indexi[0][0])
    return CijkNum

def num2cijk(num, Ord):
    #function: convert number to cijk
    #ATTENTION: here the first element is 0, not 1
    CijkAll = gen_cijall(Ord)
    Cijk = CijkAll[num - 1]
    return Cijk

def CoefForSingleMode(CrystalType, Ord, Strain):
    #function: generate the coefficient for a specified strain and crystal type
    #
    #Cijk = print_cijk(CrystalType, Ord)
    #print Cijk
    m = Strain.shape
    if len(m) == 1:
        m = 1
    else:
        m = m[0]
    StrainMode = np.zeros((m, 3, 3))
    StrainModeCoef = CoefStr(Ord)
    for i in range(2, int(Ord) + 1):
        CijUniq = np.array(symdata.get_unique(symdata.coef_ind(CrystalType, i)))
        n_uniq = len(CijUniq)
        exec('StrainModeCoef.coef' + str(i) + ' = np.zeros((m, n_uniq), dtype = float)')
    #StrainModeCoef = CoefStr(Ord)
    Cijk = CoefStr(Ord)
    for i in range(2, int(Ord) + 1):
        exec('Cijk.coef' + str(i) + ' = print_cijk(CrystalType, i)')
        for j in range(0, m):
            StrainM = crydef.vec2matrix(Strain[j, :])
            StrainMode[j, :, :] = StrainM
            #print StrainM
            exec('StrainModeCoef.coef' + str(i) + '[j, :] = gen_cijk_coef(CrystalType, StrainM, i)')
        exec('print(StrainModeCoef.coef' + str(i) + ')')
    return (Cijk, StrainModeCoef, StrainMode)

def print_cijk(CrystalType, Ord):
    #function: print the Cijk according to the crystal type and order
    CijkInd = symdata.coef_ind(CrystalType, Ord)
    CijkUniq = np.array(symdata.get_unique(CijkInd))
    Cijk = num2cijk(CijkUniq, Ord)
    print(Cijk)
    return Cijk

class CoefStr:
    """This is the structure for coefficients. The attaches is coef + i, 
    where i is 2 to Ord
    Take Ord = 3 as an example, there are two attaches, coef2 and coef3"""
    def __init__(self, Ord = 3):
        for i in range(2, int(Ord) + 1):
            exec('self.coef' + str(i) + ' = []')      
'''
a = np.zeros((4, 5))
(m, n) = a.shape
print m
'''
'''
CrystalType = 'c1'
Ord = 4

(StrainModeCoef, StrainMode) = gen_strain_mode(CrystalType, Ord)
print StrainMode
print StrainModeCoef.coef3
print StrainModeCoef.coef2
print StrainModeCoef.coef4
'''