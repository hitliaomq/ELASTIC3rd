#!python
#
import numpy as np
import elastic3rd.post.post as espost
#import elastic3rd.post.post as espost

V0 = 45.3636510818
Flag_Sep = 0
Flag_Fig = 0
Flag_Ord = 5
Flag_Data = 0
Num_Str = 1

#'''
EEnergy = "EEnergy.txt"
INPUT = "INPUT.e3rd"
STRAINMODE = "STRAINMODE"
#C = espost.post_mode(V0, Flag_Sep, Flag_Fig, Flag_Ord, EEnergy, INPUT, STRAINMODE, Num_Str)
C = espost.post_mode(V0, Flag_Sep, Flag_Fig, Flag_Data, Flag_Ord, EEnergy, INPUT, 
    STRAINMODE, Num_Str)
print(C.C2)
print(C.C3)
#print C.C4
#'''

'''
data = np.zeros((5, 9))
i_count = 0
for i in range(0, 5):
    EEnergy = 'EEnergy-C' + str(i) + '.txt'
    INPUT = 'INPUT-C' + str(i)
    STRAINMODE = 'STRAINMODE-C' + str(i)
    C = espost.post_mode(V0, Flag_Sep, Flag_Fig, Flag_Ord, EEnergy, INPUT, STRAINMODE, Num_Str)
    data[i_count, :] = np.concatenate((C.C2, C.C3), axis=0)
    i_count = i_count + 1
#np.savetxt('data-HfCN-' + str(Flag_Ord) + '.txt', data)
'''

'''
N = 11
data = np.zeros((N, 9))
i_count = 0
for i in np.arange(200, 530, 30):
    print(i)
    EEnergy = str(i) + "-EEnergy.txt"
    INPUT = str(i) + "-INPUT"
    STRAINMODE = str(i) + "-STRAINMODE"
    C = espost.post_mode(V0, Flag_Sep, Flag_Fig, Flag_Ord, EEnergy, INPUT, STRAINMODE, Num_Str)
    data[i_count, :] = np.concatenate((C.C2, C.C3), axis=0)
    print C.C2
    print C.C3
    i_count = i_count + 1
np.savetxt('data-cutoff-' + str(Flag_Ord) + '.txt', data)
'''

'''
EEnergyAll = "EEnergy-Range.txt"
INPUT_All = "INPUT-Range"
STRAINMODE = "STRAINMODE-Range"

(CrystalType, Ord, flag_se, StrainList) = espost.get_post_param(INPUT_All)
N = 13
NAll = len(StrainList)
EAll = espost.read_e(EEnergyAll)
data = np.zeros((NAll - N + 1, 9))
for i in range(0, NAll - N + 1):
    ind_s = i
    ind_e = i + N
    Ei = EAll[ind_s:ind_e, :]
    EEnergyi = "EEnergy-Range" + str(i) + ".txt"
    INPUTi = "INPUT-Base"
    np.savetxt(EEnergyi, Ei)
    input_sets = {'MaxStrain' : StrainList[ind_e-1], 'MinStrain' : StrainList[ind_s]}
    cvtest.check_input(input_sets)
    C = espost.post_mode(V0, Flag_Sep, Flag_Fig, Flag_Ord, EEnergyi, INPUTi, STRAINMODE, Num_Str)
    data[i, :] = np.concatenate((C.C2, C.C3), axis=0)
np.savetxt('data-' + str(Flag_Ord) + '-Range-6.txt', data)
'''
