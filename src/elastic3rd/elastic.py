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
#      Develop the strain-stress method for TOECs, the Flag_SE is un-used

import crystal.deform as crydef
import crystal.lattice as crylat
import symmetry.symmetry as essym
import post.post as espost
import esutils
import energy.glue
import numpy as np
import shutil
import sys
import os

INPUT = "INPUT"
ParaIn = esutils.read_input(INPUT)
eglue = __import__("energy." + ParaIn['EnergyCode'], fromlist = ParaIn['EnergyCode'])

#CrystalType = "cubic"

def elastic3(INPUT = "INPUT"):
    #Print LOGO
    esutils.print_logo()

    #Read INPUT
    #ParaIn = esutils.read_input(INPUT)    
    print("===================The input parameters===================")
    esutils.print_parain(ParaIn)

    flag_se = ParaIn['FlagSE'].lower()

    #Import glue as eglue

    #eglue = __import__("energy." + ParaIn['EnergyCode'], fromlist = ParaIn['EnergyCode'])

    #Deform mode this is for cubic
    #TODO : for any symmtry
    CrystalType = ParaIn['CrystalType']
    Ord = ParaIn['Ord']
    CalMode = ParaIn['CalMode'].lower()
    #StrainMode = []
    if CalMode == 's':
        StrainIn = esutils.read_strainmode()
        if flag_se == "e":
            Cij_mode, coef_e, StrainMode = essym.CoefForSingleMode(CrystalType, Ord, StrainIn)
        elif flag_se == "s":
            pass
    elif CalMode == 'w':
        if flag_se == "e":
            coef_e, StrainMode = essym.gen_strain_mode(CrystalType, Ord)
        elif flag_se == "s":
            pass

    #Mode_index = range(1, 7)

    StrainList = esutils.gen_strain_list(ParaIn)
    n_Strain = len(StrainList)
    n_Mode = StrainMode.shape[0]
    #n_Mode = len(Mode_index)
    E = np.zeros((n_Strain, n_Mode))
    if flag_se == "s":
        E = np.zeros((n_Strain, 6*n_Mode))
    if flag_se == "e":
        Efilename_Mode = "Energy_Mode.txt"
    elif flag_se == "s":
        Efilename_Mode = "Stress_Mode.txt"

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
    dstpath = BaseName + "/Mode0"

    E0 = get_strain_e(ParaIn, dstpath, StrainMode, BaseVec)
    if CalMode == 's':
        shutil.copyfile("STRAINMODE", BaseName + "/STRAINMODE")
    #print("\n")

    print("\n==================Deformed Crystal========================")
    
    #print E
    for i in range(1, n_Mode + 1):
        print("----------------------------------------------------------")
        print("Start calculating Mode " + str(i))
        flag_continue = 0
        ModePath = BaseName + "/Mode" + str(i)
        if ParaIn['Continue']:
            flag_continue = esutils.iscontinue(ModePath, "Mode", flag_se)
        if flag_continue:
            E = get_continue_mode_e(ModePath, E, i, flag_se)
        else:
            print(os.getcwd())
            esutils.creat_folder(ModePath)
            #os.chdir(ModePath)
            for j in range(0, n_Strain):
                print("Start calculating Strain " + str(StrainList[j]) + " in Mode " + str(i))
                if j == int(n_Strain/2):
                    if flag_se == "e":
                        E[j, i-1] = E0[0]
                    elif flag_se == "s":
                        pass
                else:
                    StrainPath = ModePath + "/Strain" + str(StrainList[j])
                    Strain = StrainList[j]/100.0
                    #print os.getcwd()
                    Eij = get_strain_e(ParaIn, StrainPath, StrainMode, BaseVec, Strain, i)
                    if flag_se == "e":
                        E[j, i-1] = Eij[0]
                        print("Energy:")
                        esutils.print_e(Eij)
                    elif flag_se == "s":
                        pass
                    os.chdir("../../../")
                print("End of Strain " + str(StrainList[j]) + " in Mode " + str(i) + "\n")
            os.chdir(ModePath)
            if flag_se == "e":
                np.savetxt(Efilename_Mode, E[:, i-1])
            os.chdir("../../")
        print("End of Mode " + str(i) + "\n")
    np.savetxt(BaseName + "/EEnergy.txt", E)
    print("\n==================Post Processing========================")
    coef_fit = espost.get_coef(StrainList/100., E, V0, flag_se, 3)
    print(coef_fit)
    coef2 = coef_e.coef2
    coef3 = coef_e.coef3
    (C2, C3) = espost.get_cij(coef_fit, coef2, coef3, flag_se)
    essym.print_cijk(CrystalType, 2)
    print(C2)
    essym.print_cijk(CrystalType, 3)
    print(C3)
    print("========================!!!END!!!=========================")


def getparam(ParaIn, dstpath, StrainOrMode):
    flag_se = ParaIn['FlagSE'].lower()
    flag_continue = 0
    if ParaIn['Continue']:
        flag_continue = esutils.iscontinue(dstpath, StrainOrMode, flag_se)    
    BaseName = ParaIn['BaseName']
    PreName = ''
    if flag_se == "e":
        PreName = "Energy"
    elif flag_se == "s":
        PreName = "Stress"
    RunStr = ''
    if ParaIn['EnergyRun']:
        RunStr = eglue.run(ParaIn['NP'], BaseName)
    else:
        RunStr = energy.glue.run()
    return (flag_se, flag_continue, BaseName, RunStr, PreName)

def get_continue_mode_e(ModePath, E, Modei, flag_se = "e"):
    os.chdir(ModePath)
    Efilename = ''
    if flag_se == "e":
        Efilename = "Energy_Mode.txt"
        print("The energy of Mode" + str(Modei) + " was calculated previously.")
        print("The energy is :")
    elif flag_se == "s":
        pass
    Etmp = np.loadtxt(Efilename)
    #E[:, i-1] = np.loadtxt(Efilename)
    print(Etmp)
    os.chdir("../../")
    if flag_se == "e":
        E[:, int(Modei-1)] = Etmp
    elif flag_se == "s":
        pass
    return E

def get_continue_strain_e(Efilename, flag_se, strain, Modei):
    if flag_se == "e":
        if Modei == 0.:
            print("The energy of undistorted structure was calculated previously.")
        else:        
            print("The energy of Strain " + str(strain) + " in Mode " + str(Modei) + "was calculated previously.")
    elif flag_se == "s":
        pass
    print(Efilename)
    E0 = np.loadtxt(Efilename)
    return E0

def get_org_strain_e(BaseName, flag_se):
    E0 = 0
    if flag_se == "e":
        E0 = eglue.get_energy(BaseName)
        print("Energy:")
        esutils.print_e(E0)
    elif flag_se == "s":
        pass
    return E0

def get_strain_e(ParaIn, dstpath, StrainMode, BaseVec, strain = 0.0, Modei = 0):
    StrainOrMode = "Strain"
    if Modei == 0:
        StrainOrMode = "Mode"
    (flag_se, flag_continue, BaseName, RunStr, PreName) = getparam(ParaIn, dstpath, StrainOrMode)
    #BaseVec = eglue.get_base_vec(BaseName)
    Efilename = PreName + "_" + StrainOrMode + ".txt"

    if flag_continue:        
        os.chdir(dstpath)
        #E0 is energy for strain-energy method, and stress for strain-stress method
        E0 = get_continue_strain_e(Efilename, flag_se, strain, Modei)
        if Modei == 0:
            os.chdir("../../")
    else:
        if Modei == 0:
            esutils.creat_folder(BaseName)
            shutil.copyfile("INPUT", BaseName + "/INPUT")
        #print os.getcwd()
        esutils.creat_folder(dstpath)
        eglue.copy_files(BaseName, dstpath)
        os.chdir(dstpath)
        if Modei > 0:
            StrainMatrix = strain*StrainMode[Modei-1]
            BaseVecNew = BaseVec .dot (crydef.strain2deformgrad(StrainMatrix))
            crylat.print_lattice(BaseVecNew) 
            eglue.write_base_vec(BaseName, BaseVecNew)
        os.system(RunStr + ">> FP_OUT")
        E0 = get_org_strain_e(BaseName, flag_se)
        print(E0)
        np.savetxt(Efilename, E0)    
        if Modei == 0:
            os.chdir("../../")
    return E0

# end 
