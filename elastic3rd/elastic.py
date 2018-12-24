#!python
#the workflow of elastic3rd
# 1. calculate the energy of undeformed structure
# 2. deform the crystal by deform mode
# 3. calculate the energy of deformed structure
# 4. fit the function of energy vs strain, get the coefficient of fitting
# 5. solve the corresponding equation, get the 2nd and 3rd elastic stiffness
#
#The result tree
# BaseName
# 
#TODO: find the symmtry and determine the deform mode

import crystal.deform as crydef
import crystal.lattice as crylat
import post.post as espost
import esutils
import energy.glue
import numpy as np
import sys
import os

#Print LOGO
esutils.print_logo()

#Read INPUT
if len(sys.argv) > 1:
    INPUT = sys.argv[1]
else:
    INPUT = "INPUT"
ParaIn = esutils.read_input(INPUT)
crystaltype = "cubic"
print("===================The input parameters===================")
esutils.print_parain(ParaIn)

#Import glue as eglue

eglue = __import__("energy." + ParaIn['EnergyCode'], fromlist = ParaIn['EnergyCode'])

def elastic3():
    #Deform mode this is for cubic
    #TODO : for any symmtry
    Mode_index = range(1, 7)

    #Get the energy run string
    BaseName = ParaIn['BaseName']
    if ParaIn['EnergyRun']:
        RunStr = eglue.run(ParaIn['NP'], BaseName)
    else:
        RunStr = energy.glue.run()

    #os.system(RunStr)
    #Get the base vector
    BaseVec = eglue.get_base_vec(BaseName)
    print("====================Crystal Structure====================")
    V0 = crylat.print_lattice(BaseVec)

    #Calculate the energy of undeformed structure
    #    creat the folder Mode0
    esutils.creat_folder(BaseName)
    dstpath = BaseName + "/Mode0"
    esutils.creat_folder(dstpath)
    eglue.copy_files(BaseName, dstpath)
    os.chdir(dstpath)
    os.system(RunStr)
    E0 = eglue.get_energy(BaseName)
    os.chdir("../../")
    print("Energy of undeformed structure:")
    esutils.print_e(E0)
    #print("\n")

    print("\n==================Deformed Crystal========================")
    StrainList = esutils.gen_strain_list(ParaIn)
    n_Strain = len(StrainList)
    n_Mode = len(Mode_index)
    E = np.zeros((n_Strain, n_Mode))
    print E
    for i in range(1, n_Mode + 1):
        print("----------------------------------------------------------")
        print("Start calculating Mode " + str(i))
        ModePath = BaseName + "/Mode" + str(i)
        esutils.creat_folder(ModePath)
        for j in range(0, n_Strain):
            print("Start calculating Strain " + str(StrainList[j]) + " in Mode " + str(i))
            if j == int(n_Strain/2):
                E[j][i-1] = E0[0]
            else:
                StrainPath = ModePath + "/Strain" + str(StrainList[j])
                esutils.creat_folder(StrainPath)
                eglue.copy_files(BaseName, StrainPath)
                os.chdir(StrainPath)
                Strain = StrainList[j]/100.0
                StrainVerctor = crydef.deform_mode(Strain, i)
                StrainMatrix = crydef.vec2matrix(StrainVerctor)
                BaseVecNew = BaseVec .dot (crydef.strain2deformgrad(StrainMatrix))
                eglue.write_base_vec(BaseName, BaseVecNew)
                os.system(RunStr)
                Eij = eglue.get_energy(BaseName)                
                E[j][i-1] = Eij[0]
                crylat.print_lattice(BaseVecNew)
                print("Energy:")
                esutils.print_e(Eij)
                os.chdir("../../../")
            print("End of Strain " + str(StrainList[j]) + " in Mode " + str(i) + "\n")
        print("End of Mode " + str(i) + "\n")
    np.savetxt(BaseName + "/EEnergy.txt", E)
    print("\n==================Post Processing========================")
    coef_fit = espost.get_coef(StrainList, E, V0)
    print coef_fit
    (C2, C3) = espost.get_cij(coef_fit, crystaltype)
    print C2
    print C3
    print("========================!!!END!!!=========================")

# end 
