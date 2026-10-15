import numpy as np
import elastic3rd.esutils as esutils
import elastic3rd.post.post as espost
import elastic3rd.symmetry.symmetry as essym
import elastic3rd.symmetry.symmetry_stress as symstr
import elastic3rd.post.esfit as esfit
import elastic3rd.tools.comtools as comtools

'''
V0 = 45.3636510818
Flag_Sep = 1
Flag_Fig = 0
Flag_Ord = 7
Num_Str = 1
Nsingle = 21   # The number of points used in a single fitting

EEnergy = "EEnergy.txt"
INPUT = "INPUT"
STRAINMODE = "STRAINMODE"
'''

def post_asym(INPUT, STRAINMODE, EEnergy, Nsingle, Flag_Ord = 3, Num_Str = 1, V0 = 0, Flag_Fig = 0, Flag_Sep = 0):
    E = np.loadtxt(EEnergy)
    StrainIn = esutils.read_strainmode(STRAINMODE)
    (CrystalType, Ord, flag_se, StrainList) = espost.get_post_param(INPUT)
    ParaIn = esutils.read_input(INPUT)
    BaseName = ParaIn['BaseName']
    resultFolder = BaseName + "-asym-Result"
    esutils.creat_folder(resultFolder)
    C = espost.CCOEF(Ord)
    Nall = len(StrainList)

    C2 = np.zeros([Nall + 1 - Nsingle, 4])
    C3 = np.zeros([Nall + 1 - Nsingle, 7])
    for i in range(Nsingle, Nall + 1):
        Ei = E[i - Nsingle:i, :]
        Straini = StrainList[i - Nsingle:i]
        C2[i - Nsingle, 0] = (Straini[0] + Straini[-1])/2.
        C3[i - Nsingle, 0] = (Straini[0] + Straini[-1])/2.
        if flag_se == "e":
            Cij_mode, coef_e, StrainMode = essym.CoefForSingleMode(CrystalType, Ord, StrainIn)
        elif flag_se == "s":
            Cij_mode, coef_e, StrainMode = symstr.CoefForSingleMode(CrystalType, Ord, StrainIn)
        coef_fit = espost.get_coef(Straini/100., Ei, V0, flag_se, Flag_Ord)
        if Flag_Fig == 1:
            esfit.multiesplot(Straini/100., Ei, coef_fit, flag_se, Flag_Ord, V0)
        n = eval('coef_e.coef' + str(int(Ord)) + '.shape[1]')
        SMRank = eval('np.linalg.matrix_rank(coef_e.coef' + str(int(Ord)) + ')')
        if SMRank < n:
            print(SMRank)
            print(n)
            print('The strain mode is not full-rank.')
        else:
            C = espost.get_cijall(coef_fit, coef_e, int(Ord), flag_se)
            C2[i - Nsingle, 1:] = C.C2
            C3[i - Nsingle, 1:] = C.C3       
            
    np.savetxt(resultFolder + "/" + BaseName + "-asym-C2-Ord" + str(Flag_Ord) + ".txt", C2)
    np.savetxt(resultFolder + "/" + BaseName + "-asym-C3-Ord" + str(Flag_Ord) + ".txt", C3)
    return BaseName

def postn(INPUT, STRAINMODE, EEnergy, Nsingle, n = 7, Num_Str = 1, V0 = 0, Flag_Fig = 0, Flag_Sep = 0):
    for i in range(3, n+1):
        BaseName = post_asym(INPUT, STRAINMODE, EEnergy, Nsingle, i, Num_Str, V0, Flag_Fig, Flag_Sep)
    flag = "asym"
    comtools.sortByOrder(BaseName, flag, n, OrdMin = 3)