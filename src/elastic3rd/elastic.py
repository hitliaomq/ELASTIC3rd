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
import importlib

#INPUT = "INPUT"
#ParaIn = esutils.read_input(INPUT)
#eglue = __import__("energy." + ParaIn['EnergyCode'], fromlist = ParaIn['EnergyCode'])

#CrystalType = "cubic"

def elastic3(INPUT = "INPUT"):
    #This is the main function for calculation
    
    #Print LOGO
    esutils.print_logo()

    #Read INPUT
    ParaIn = esutils.read_input(INPUT)    
    print("===================The input parameters===================")
    esutils.print_parain(ParaIn)

    flag_se = ParaIn['FlagSE'].lower()
    EnergyCode = ParaIn['EnergyCode'].lower()

    #Import first principles code module according to the EnergyCode settings in INPUT file
    eglue = importlib.import_module('elastic3rd.energy.{}'.format(EnergyCode))

    #Get some parameter's value from INPUT
    CrystalType = ParaIn['CrystalType']
    Ord = ParaIn['Ord']
    CalMode = ParaIn['CalMode'].lower()
    
    #Get the strain modes according to symmetry and order, the strain-stress method is under development
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

    #Get the strain list
    StrainList = esutils.gen_strain_list(ParaIn)
    n_Strain = len(StrainList)
    n_Mode = StrainMode.shape[0]

    #Init the energy
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

    #Begin the calculation of deformed structure
    print("\n==================Deformed Crystal========================")

    #Do a loop over the strain modes
    for i in range(1, n_Mode + 1):
        print("----------------------------------------------------------")
        print("Start calculating Mode " + str(i))
        flag_continue = 0

        #Judge if this mode has been calculated or not
        ModePath = BaseName + "/Mode" + str(i)
        if ParaIn['Continue']:
            flag_continue = esutils.iscontinue(ModePath, "Mode", flag_se)
        if flag_continue:
            #Calculated previously
            E = get_continue_mode_e(ModePath, E, i, flag_se)
        else:
            #Not calculated yet, new calculations
            print(os.getcwd())
            esutils.creat_folder(ModePath)
            #Do a loop over strain list
            for j in range(0, n_Strain):
                print("Start calculating Strain " + str(StrainList[j]) + " in Mode " + str(i))
                if j == int(n_Strain/2):
                    #This is the un-deformed structure
                    if flag_se == "e":
                        E[j, i-1] = E0[0]
                    elif flag_se == "s":
                        pass
                else:
                    #Deformed structure
                    StrainPath = ModePath + "/Strain" + str(StrainList[j])
                    Strain = StrainList[j]/100.0
                    #The judge on the continue mode is done in the following function
                    Eij = get_strain_e(ParaIn, StrainPath, StrainMode, BaseVec, Strain, i)
                    if flag_se == "e":
                        E[j, i-1] = Eij[0]
                        print("Energy:")
                    elif flag_se == "s":
                        pass
                    os.chdir("../../../")
                print("End of Strain " + str(StrainList[j]) + " in Mode " + str(i) + "\n")
            os.chdir(ModePath)
            if flag_se == "e":
                np.savetxt(Efilename_Mode, E[:, i-1])
            os.chdir("../../")
        print("End of Mode " + str(i) + "\n")

    #Sve the energy
    np.savetxt(BaseName + "/EEnergy.txt", E)

    ##Post process
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
    '''
    Get the parameters for calculations
    Parameters
    ----------
        ParaIn: dict
            The input parameter dict
        dstpath: str
            The destination path
        StrainOrMode: str
            Strain or Mode
    Return
    ------
        flag_se: str
            e or s. e for strain-energy method, s for strain-stress method
        flag_continue: 0 or 1
            0 or 1. 0 for current mode or strain is not calculated yet, 1 for calculated
        BaseName: str
            The basename, is taken from ParaIn directely
        RunStr: str
            Get the string for running the first principles code.
        PreName: str
            The pre of the filename of result file. It is contraled by the flag_se.
                If flag_se='e', PreName='Energy', else PreName='Stress'
    '''
    flag_se = ParaIn['FlagSE'].lower()
    EnergyCode = ParaIn['EnergyCode'].lower()

    eglue = importlib.import_module('elastic3rd.energy.{}'.format(EnergyCode))
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
    '''
    Get the energy or stress of previous mode folder
    Parameters
    ----------
        ModePath: str
            The path or the current
        E: np.ndarray, size: n_strain x n_mode
            The initial of result. It can be determined by the number of strainmodes and steps.
        Modei: int
            The ith mode
        flag_se: str
            s or e. e for strain-energy method and s for strain-stress method
    Return
    ------
        E: np.ndarray
            The updated E
    '''
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
    '''
    Get the energy or stress of previous mode folder
    Parameters
    ----------
        Efilename: str
            The filename of the energy file
        strain: float
            the strain, e.g. 3.0
        Modei: int
            The ith mode
        flag_se: str
            s or e. e for strain-energy method and s for strain-stress method
    Return
    ------
        E0: np.ndarray
            The energy of current strain
    '''
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

def get_org_strain_e(BaseName, flag_se='e', EnergyCode='vasp'):
    '''
    Get the energy from first principles codes and pring it.
    Parameters
    ----------
        BaseName: str
            The basename, the name of the cell file or param file(for CASTEP)
        flag_se: str
            s or e. e for strain-energy method and s for strain-stress method
    Return
    ------
        E0: list
            The energy of current strain. in multi unit
    '''
    eglue = importlib.import_module('elastic3rd.energy.{}'.format(EnergyCode))
    E0 = 0
    if flag_se == "e":
        E0 = eglue.get_energy(BaseName)
        print("Energy:")
        esutils.print_e(E0)
    elif flag_se == "s":
        pass
    return E0

def get_strain_e(ParaIn, dstpath, StrainMode, BaseVec, strain=0.0, Modei=0):
    '''
    In the strain folder, run first principles code and get the energy
    Parameters
    ----------
        ParaIn: dict
            The dict of input parameters
        dstpath: str
            The destination path 
        StrainMode: 3D np.ndarray, size: nx3x3
            The strain mode.
        BaseVec: np.ndarray
            The crystal vector
        strain: float
            The value of the strain. e.g. 0.03
        Modei: int
            The ith strain mode. When modei=0, it means the undeformed structure
    Return
    ------
        E0: list
            The energy of current strain. in multi unit
    '''
    StrainOrMode = "Strain"
    if Modei == 0:
        StrainOrMode = "Mode"
    (flag_se, flag_continue, BaseName, RunStr, PreName) = getparam(ParaIn, dstpath, StrainOrMode)
    EnergyCode = ParaIn['EnergyCode'].lower()

    eglue = importlib.import_module('elastic3rd.energy.{}'.format(EnergyCode))
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
        E0 = get_org_strain_e(BaseName, flag_se, EnergyCode=EnergyCode)
        np.savetxt(Efilename, E0)    
        if Modei == 0:
            os.chdir("../../")
    return E0

# end 
