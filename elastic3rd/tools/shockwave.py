#!python
#
#This method is proposed in PRB 93, 174113 (2016)
#  In this method, it is comparable with experiments 
#  However, the higher order effect is not taken into consideration
#  Hence, there are some errors
#  Ref. [Thesis] Lang, John Michael. 
#         Mechanical and Optical Response of Diamond Crystals Shock Compressed along Different Orientations
#         Washington State University

'''
direction  a    b    c      STRAINMODE
100       100   010  001    1 0 0 0 0 0   P.18 in Thesis
110       110   -110 001    0.5 0.5 0 0 0 1   P.19 in Thesis
111       111   1-10 11-2   0.3 0.3 0.3 0.6 0.6 0.6   P.20 in Theis
The STRAINMODE is stored in STRAINMODE-COMPRESSION
The strain amplitude is always negtive, corresponding to compression

Transformation matrix
100        [[1., 0, 0], [0, 1., 0], [0, 0, 1.]]
110        [[1./sqrt(2.), 1./sqrt(2.), 0], [-1./sqrt(2.), 1./sqrt(2.), 0], [0, 0, 1.]]
111        [[1./sqrt(3.), 1./sqrt(3.), 1./sqrt(3.)], [-1./sqrt(6.), -1./sqrt(6.), 2./sqrt(6.)], [1./sqrt(2.), -1./sqrt(2.), 0]]
'''
import os
import math
import numpy as np
import elastic3rd
import elastic3rd.esutils as esutils
import elastic3rd.crystal.deform as crydef
import elastic3rd.crystal.lattice as crylat
import elastic3rd.symmetry

flag_test = 0

INPUT = "INPUT.e3rd"
if os.path.exists(INPUT):
    ParaIn = esutils.read_input(INPUT)
    eglue = __import__("elastic3rd.energy." + ParaIn['EnergyCode'], fromlist = ParaIn['EnergyCode'])
    print(ParaIn['EnergyCode'])
else:
    eglue = __import__("elastic3rd.energy.vasp", fromlist = "vasp")
    print("vasp")

BaseName = ParaIn['BaseName']
BaseVec = eglue.get_base_vec(BaseName)
V0 = crylat.print_lattice(BaseVec)
StrainList = esutils.gen_strain_list(ParaIn)
sqrt2 = math.sqrt(2.)
sqrt3 = math.sqrt(3.)
sqrt6 = math.sqrt(6.)
TransfMat = {"100" : [[1., 0, 0], [0, 1., 0], [0, 0, 1.]],
             "110" : [[1./sqrt2, 1./sqrt2, 0], [-1./sqrt2, 1./sqrt2, 0], [0, 0, 1.]],
             "111" : [[1./sqrt3, 1./sqrt3, 1./sqrt3], [-1./sqrt6, -1./sqrt6, 2./sqrt6], [1./sqrt2, -1./sqrt2, 0]]}
NewAxesall = {"100" : np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 
              "110" : np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]]),
              "111" : np.array([[1, 1, 1], [1, -1, 0], [1, 1, -2]])}
StrainMode = {"100" : np.array([1., 0, 0, 0, 0, 0]), 
              "110" : np.array([0.5, 0.5, 0, 0, 0, 1.]), 
              "111" : np.array([0.3, 0.3, 0.3, 0.6, 0.6, 0.6])}
SMOrder = ["100", "110", "111"]

def post_shockwave(SMOrder, BaseName, StrainList):
    
    for i in range(0, len(SMOrder)):
        #Read the stress of initial structure
        Direct = SMOrder[i]
        StrainModei = StrainMode[Direct]
        FilePath = BaseName + "/Mode0"
        cauchy0 = readAndTrans(FilePath, Direct, Mode = 1)
        FolderName = BaseName + "/Mode" + str(int(i+1))
        OutFileName = BaseName + "-shockwave-" + str(i) + ".txt"
        Output = np.zeros((len(StrainList), 7))
        i_count = 0
        for j in StrainList:
            straini = j/100.
            StrainPrem = getStrainFromRho(straini, StrainModei, BaseVec)
            if int(j) == 0:
                cauchy = cauchy0
            else:
                FilePath = FolderName + "/Strain" + str(j)
                cauchy = readAndTrans(FilePath, Direct, Mode = 0)
            Output[i_count, 0] = StrainPrem
            V = getV(straini, StrainModei, BaseVec)
            Output[i_count, 1:-1] = (V0/V) * stresstensor2vec(cauchy)
            i_count = i_count + 1
        np.savetxt(OutFileName, Output)

def stresstensor2vec(S):
    Stress = np.array([S[0,0], S[1,1], S[2,2], S[1, 2], S[0, 2], S[0, 1]])
    return Stress

def getStrainFromRho(strain, StrainMode, BaseVec):
    V = getV(strain, StrainMode, BaseVec)
    StrainPrem = 0.5*((V/V0)**2 - 1)
    return StrainPrem

def getV(strain, StrainMode, BaseVec):
    StrainMatrix = strain*StrainMode
    BaseVecNew = BaseVec .dot (crydef.strain2deformgrad(StrainMatrix))
    V = crylat.print_lattice(BaseVecNew)
    return V

def readAndTrans(FilePath, Direct, Mode = 0):
    if Mode == 1:
        FileName = FilePath + "/Stress_Mode.txt"
    else:
        FileName = FilePath + "/Stress_Stress.txt"
    Stress = np.loadtxt(FileName)
    NewAxesi = NewAxesall[Direct]
    StressP = StressTransfer(Stress, NewAxesi)
    return StressP

def run_strainmode_single():
    pass

def StressTransfer(Stress, NewAxes, OriAxes = np.array([[1., 0, 0], [0, 1., 0], [0, 0, 1.]])):
    TMatrix = np.einsum("ij, jk", OriAxes, NewAxes)
    # Normalize TMatrix
    print(TMatrix)
    TMatrix = TMatrix / np.reshape(np.sqrt(np.einsum("ij->i", np.einsum("ij,ij->ij", TMatrix, TMatrix))), (-1, 1))
    print(TMatrix)
    StressP = np.einsum("ik, jl, kl->ij", TMatrix, TMatrix, Stress)
    return StressP

def CalcTOECs(C2nd, dC2dP, Cprem):
    # PHYSICAL REVIEW B 93, 174113 (2016)
    a10 = -dC2dP[0]*(C2nd[0] + 2.*C2nd[1])-C2nd[0]
    C111 = Cprem[0]
    C112 = (a10 - C111)/2.
    a11 = -dC2dP[1]*(C2nd[0] + 2.*C2nd[1])-C2nd[1]
    C123 = (a11 - 2.*C112)
    #C123 = (a11 - C112)/2.  #prb-2016
    C166 = (4.*(Cprem[1])-C111-3.*C112)/12
    a12 = -dC2dP[2]*(C2nd[0] + 2.*C2nd[1])-C2nd[2]
    C144 = a12 - 2.*C166
    C456 = (9.*Cprem[2] - C111 - 6.*C112 - 2.*C123 - 12.*C144 - 24.*C166)/16
    TOECs = np.array([C111, C112, C123, C144, C166, C456])
    return TOECs

'''
dC2dP = np.array([10.61724, -0.40227, 1.26163])
C2nd = np.array([1054, 124, 559])
Cprem = np.array([-7515, -15183, -14385])
TOECs = CalcTOECs(C2nd, dC2dP, Cprem)
print(TOECs)
'''