#!python
#the workflow of elastic3rd
# 1. Geom optimize the initial structure
# 2. Add a rand offset or a small perturbation
# 3. Calculate the energy of initial structure
# 4. Deform the crystal by deform mode
# 5. Calculate the energy of deformed structure
# 6. Fit the function of energy vs strain, get the coefficient of fitting
# 7. Solve the corresponding equation, get the 2nd and 3rd elastic stiffness
#
#The result tree
# BaseName
# 
#TODO: find the symmtry and determine the deform mode
#author: Mingqing Liao
#e-mail: liaomq1900127@163.com
#FGMS group @ Harbin Institute of Technology(HIT)
#Mingqing Liao @ Jiangsu University of Science and Technology


import logging
import elastic3rd.crystal.deform as crydef
import elastic3rd.crystal.lattice as crylat
import elastic3rd.symmetry.symmetry as essym
import elastic3rd.symmetry.symmetry_stress as symstr
import elastic3rd.symmetry.gencoef as gencoef
import elastic3rd.post.post as espost
import elastic3rd.post.post2 as espost2
import elastic3rd.post.post_nd as post_nd
import elastic3rd.esutils as esutils
from elastic3rd.esutils import logger
import elastic3rd.energy.glue
import numpy as np
import shutil
import sys
import os
import importlib

RUNDIRECTCODE = ['mattersim', 'chgnet', 'matgl', 'sevenn', 
                 'orb', 'mace', 'grace', 'nequip', 'dpmd', 'alphanet', 'tace', 'matris']

def import_energycode(EnergyCode="vasp", INPUT="INPUT.e3rd"):
    if os.path.exists(INPUT):
        ParaIn = esutils.read_input(INPUT)
        eglue = importlib.import_module('energy.{}'.format(ParaIn['EnergyCode']))
        #eglue = __import__("energy." + ParaIn['EnergyCode'], fromlist = ParaIn['EnergyCode'])
        logger.info(ParaIn['EnergyCode'])
    else:
        eglue = importlib.import_module('energy.vasp')
        #eglue = __import__("energy.vasp", fromlist = "vasp")
        logger.info("vasp")
    return eglue

def elastic3(INPUT = "INPUT.e3rd"):
    #logger = init_log()
    #esutils.GetNowTime()
    #Print LOGO
    esutils.print_logo()

    ################################# READ INPUT #########################################
    ParaIn = esutils.read_input(INPUT)
    GeomMap = {1: ['G'], 2: ['GS', 'G'], 0: 0, 31: ['G'], 32: ['GS', 'G']}
    logger.info("===================The input parameters===================")
    esutils.print_parain(ParaIn)

    flag_se = ParaIn['FlagSE'].lower()
    coef_method = ParaIn['CoefMethod'].lower()
    acc = ParaIn['Acc']
    CrystalType = ParaIn['CrystalType']
    Ord = ParaIn['Ord']
    CalMode = ParaIn['CalMode'].lower()
    BaseName = ParaIn['BaseName']
    geom = ParaIn['Geom']
    Geom = GeomMap[ParaIn['Geom']]
    AddPert = ParaIn['AddPert']
    PertMethod = ParaIn['PertMethod']
    EnergyCode = ParaIn['EnergyCode'].lower()

    #load the energy module
    eglue = importlib.import_module('elastic3rd.energy.{}'.format(EnergyCode))

    #Get the strainmode from build-in database
    if (CalMode == 'w') and (ParaIn['GenSM'] == 1):
        # If use GenSM=1 and CalMode='w', generate the strainmode, then is the same as CalMode='s'
        sm_bldin = gencoef.get_strainmode(CrystalType, Ord, flag_se)
        np.savetxt("STRAINMODE", sm_bldin)
        CalMode = 's'
    
    #Generate the strain mode and corresponding coefficients
    if CalMode == 's':
        StrainIn = esutils.read_strainmode()
        if flag_se == "e":
            Cij_mode, coef_e, StrainMode = essym.CoefForSingleMode(CrystalType, Ord, StrainIn)
        elif flag_se == "s":
            Cij_mode, coef_e, StrainMode = symstr.CoefForSingleMode(CrystalType, Ord, StrainIn)
    elif CalMode == 'w':
        if flag_se == "e":
            coef_e, StrainMode = essym.gen_strain_mode(CrystalType, Ord)
        elif flag_se == "s":
            coef_e, StrainMode = symstr.gen_strain_mode(CrystalType, Ord)
    if flag_se == 's':
        CalMethod = 'Strain-Stress Method'
    else:
        CalMethod = 'Strain-Energy Method'
    logger.info("For crystal type({}), {} number of strain mode is required for {} method.".format(CrystalType, len(StrainMode), CalMethod))

    #Get the strain list
    StrainList = esutils.gen_strain_list(ParaIn)
    n_Strain = len(StrainList)
    n_Mode = StrainMode.shape[0]
    
    #Initial the energy or stress
    E = np.zeros((n_Strain, n_Mode))
    if flag_se == "s":
        E = np.zeros((n_Strain, 6*n_Mode))
    if flag_se == "e":
        Efilename_Mode = "Energy_Mode.txt"
    elif flag_se == "s":
        Efilename_Mode = "Stress_Mode.txt"

    #Get the energy run string
    if ParaIn['EnergyCode'].lower() in RUNDIRECTCODE:
        RunStr = ''
    else:
        if ParaIn['EnergyRun']:
            RunStr = eglue.run(ParaIn['NP'], BaseName)
        else:
            RunStr = elastic3rd.energy.glue.run()

    ################################## BEGIN THE CALCULATION ####################################
    #os.system(RunStr)
    #Get the base vector
    BaseVec = eglue.get_base_vec(BaseName)
    logger.info("====================Crystal Structure====================")
    esutils.GetNowTime()
    V0 = crylat.print_lattice(BaseVec)

    #Geom optimize the initial structure
    if Geom:
        for i, Geomi in enumerate(Geom):
            if Geomi == 'G':
                geom_path = "Geom"
            else:
                geom_path = "Geom-{}".format(Geomi)
            E_geom = get_strain_e(ParaIn, geom_path, StrainMode, BaseVec, Modei=-1, Geomi=Geomi, INPUT=INPUT)
            logger.info("The energy/stress of Geom is \n{}".format(E_geom))
            #Update the baseVec
            BaseVec = eglue.get_base_vec(BaseName)
        logger.info("====================Crystal Structure after Geom====================")
        V0 = crylat.print_lattice(BaseVec)
        if geom > 30:
            return

    #Add perturbation
    if AddPert:
        eglue.add_perturbation(BaseName, tor=2.0E-10, PertMethod=PertMethod)

    #Calculate the energy/stress of undeformed structure
    dstpath = os.path.join(BaseName, 'Mode0')
    E0 = get_strain_e(ParaIn, dstpath, StrainMode, BaseVec, Modei=0, INPUT=INPUT)

    if CalMode == 's':
        shutil.copyfile("STRAINMODE", os.path.join(BaseName, 'STRAINMODE'))
    if os.path.exists('STRAINAMP'):
        shutil.copyfile("STRAINAMP", os.path.join(BaseName, 'STRAINAMP'))
    
    #Begin to calculate the stress/energy of deformed structure
    esutils.GetNowTime()
    logger.info("\n==================Deformed Crystal========================\n")
    
    for i in range(1, n_Mode + 1):
        logger.info("----------------------- Mode {} ---------------------------\n".format(i))
        esutils.GetNowTime()
        flag_continue = 0
        ModePath = os.path.join(BaseName, 'Mode{}'.format(str(i)))

        #If current mode is calculated previousely
        if ParaIn['Continue']:
            flag_continue = esutils.iscontinue(ModePath, "Mode", flag_se)
        if flag_continue:
            E = get_continue_mode_e(ModePath, E, i, flag_se)
        else:
            logger.info(os.getcwd())
            esutils.creat_folder(ModePath)
            #os.chdir(ModePath)
            for j in range(0, n_Strain):
                logger.info("------------------ Strain {} in Mode {} ------------------\n".format(StrainList[j], i))
                esutils.GetNowTime()
                if StrainList[j] == 0:
                    if flag_se == "e":
                        E[j, i-1] = E0[0]
                    elif flag_se == "s":
                        #(PK, Etmp) = symstr.cauchy2pk(E0, StrainMode[i-1])
                        E[j, 6*(i-1):6*i] = np.array([E0[0,0], E0[1,1], E0[2,2], E0[1,2], E0[0, 2], E0[0,1]])
                else:
                    StrainPath = os.path.join(ModePath, 'Strain{}'.format(str(StrainList[j])))
                    Strain = StrainList[j]/100.0
                    #print os.getcwd()
                    Eij = get_strain_e(ParaIn, StrainPath, StrainMode, BaseVec, strain=Strain, Modei=i, INPUT=INPUT)
                    if flag_se == "e":
                        E[j, i-1] = Eij[0]
                        logger.info("Energy:")
                        esutils.print_e(Eij)
                    elif flag_se == "s":
                        logger.info(StrainMode[i-1])
                        (PK, Etmp) = symstr.cauchy2pk(Eij, Strain*StrainMode[i-1])
                        E[j, 6*(i-1):6*i] = Etmp
                        logger.info("Cauchy stress:")
                        logger.info(Eij)
                        logger.info("Lagrangian stress:")
                        logger.info(Etmp)                    
                    os.chdir("../../../")
                logger.info("------------------ End of Strain {} in Mode {} ------------------\n".format(StrainList[j], i))
            
            #Save the stress/energy of mode i
            os.chdir(ModePath)
            if flag_se == "e":
                np.savetxt(Efilename_Mode, E[:, i-1])
            elif flag_se == "s":
                np.savetxt(Efilename_Mode, E[:, 6*(i-1):6*i])
            os.chdir("../../")
        logger.info("------------------------ End of Mode {} --------------------------\n".format(i))
    
    #Save stress/energy of all strain modes
    np.savetxt(os.path.join(BaseName, 'EEnergy.txt'), E)

    logger.info("\n=============================== Post Processing ==================================\n")
    if flag_se == 's':
        (PK0, P0) = symstr.cauchy2pk(E0, np.zeros((3, 3)))
    else:
        P0 = np.zeros((1, 6))
    if coef_method == 'num':
        coef_fit = post_nd.get_coef_nd(StrainList/100., E, V0, flag_se, P0, acc=acc)
    else:
        coef_fit = espost.get_coef(StrainList/100., E, V0, flag_se, 3)
    logger.info(coef_fit)
    #coef2 = coef_e.coef2
    #coef3 = coef_e.coef3
    #(C2, C3) = espost.get_cij(coef_fit, coef2, coef3, flag_se)
    C = espost2.get_cij(coef_fit, coef_e, flag_se, Ord=Ord)
    essym.print_cijk(CrystalType, 2)
    logger.info(C[-1])
    if Ord > 2:
        essym.print_cijk(CrystalType, 3)
        logger.info(C[-2])
    logger.info("========================!!!END!!!=========================")
    return C

##################################### END OF MAIN ##############################################

def getparam(ParaIn, dstpath, StrainOrMode):
    '''
    Get the parameters for calculations
    '''
    eglue = importlib.import_module('elastic3rd.energy.{}'.format(ParaIn['EnergyCode'].lower()))
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
    if ParaIn['EnergyCode'].lower() in RUNDIRECTCODE:
        pass
    else:
        if ParaIn['EnergyRun']:
            RunStr = eglue.run(ParaIn['NP'], BaseName)
        else:
            RunStr = elastic3rd.energy.glue.run()
    return (flag_se, flag_continue, BaseName, RunStr, PreName)

def get_continue_mode_e(ModePath, E, Modei, flag_se="e"):
    '''
    Get the stress/energy of previous calculated result (in mode level)
    '''
    os.chdir(ModePath)
    Efilename = ''
    if flag_se == "e":
        Efilename = "Energy_Mode.txt"
        logger.info("The energy of Mode" + str(Modei) + " was calculated previously.")
        logger.info("The energy is :")
    elif flag_se == "s":
        Efilename = "Stress_Mode.txt"
        logger.info("The stress of Mode" + str(Modei) + " was calculated previously.")
        logger.info("The stress is :")
    Etmp = np.loadtxt(Efilename)
    #E[:, i-1] = np.loadtxt(Efilename)
    logger.info(Etmp)
    os.chdir("../../")
    if flag_se == "e":
        E[:, int(Modei-1)] = Etmp
    elif flag_se == "s":
        E[:, int(6*(Modei-1)):int(6*Modei)] = Etmp
    return E

def get_continue_strain_e(Efilename, flag_se, strain, Modei):
    '''
    Get the stress/energy of previous calculated result (in strain level)
    '''
    if flag_se == 'e':
        prop = 'energy'
    else:
        prop = 'stress'
    if int(Modei) == 0:
        logger.info("The {} of undistorted structure was calculated previously.".format(prop))
    else:        
        logger.info("The {} of Strain {} in Mode {} was calculated previously.".format(prop, strain, Modei))       
    E0 = np.loadtxt(Efilename)
    return E0

def get_org_strain_e(BaseName, flag_se, EnergyCode='vasp'):
    '''
    Get the stress/energy from FP result
    '''
    eglue = importlib.import_module('elastic3rd.energy.{}'.format(EnergyCode.lower()))
    E0 = 0
    if flag_se == "e":
        E0 = eglue.get_energy(BaseName)
        logger.info("Energy:")
        esutils.print_e(E0)
    elif flag_se == "s":
        E0 = eglue.get_stress(BaseName)
        logger.info("Cauchy stress (in GPa):")
        logger.info(E0)
    return E0

def get_strain_e(ParaIn, dstpath, StrainMode, BaseVec, strain=0.0, Modei=0, RemLargFile=1, Geomi='G', INPUT='INPUT.e3rd'):
    '''
    Run FP code and get the stress or energy (Geom/Single/DeformStructure)

    Parameter
    ---------
        ParaIn: dict
            The input parameter
        dstpath: str(path-like)
            The destination path of the calculation
        StrainMode: np.ndarray (nx3x3)
            The strain mode
        BaseVec: np.ndarray (3x3)
            The lattice vector
        strain: float
            The strain
        Modei: int
            -1: full geom optimization
             0: ion-relax (initial structure)
             N: ion-relax (the Nth StrainMode at strain)
    Return
    ------
    '''
    eglue = importlib.import_module('elastic3rd.energy.{}'.format(ParaIn['EnergyCode'].lower()))
    StrainOrMode = "Strain"
    if Modei == 0:
        StrainOrMode = "Mode"
    elif Modei == -1:
        StrainOrMode = "Geom"
    (flag_se, flag_continue, BaseName, RunStr, PreName) = getparam(ParaIn, dstpath, StrainOrMode)
    #BaseVec = eglue.get_base_vec(BaseName)
    Efilename = PreName + "_" + StrainOrMode + ".txt"

    if flag_continue:        
        os.chdir(dstpath)
        #E0 is energy for strain-energy method, and stress for strain-stress method
        E0 = get_continue_strain_e(Efilename, flag_se, strain, Modei)
        if Modei == 0:
            os.chdir("../../")
        elif Modei == -1:
            os.chdir("../")
    else:
        if Modei == 0:
            esutils.creat_folder(BaseName)
            shutil.copyfile(INPUT, os.path.join(BaseName, INPUT))
        #elif Modei == -1:
        #    esutils.creat_folder('Geom')
        #print os.getcwd()
        esutils.creat_folder(dstpath)
        eglue.copy_files(BaseName, dstpath)
        os.chdir(dstpath)
        if Modei > 0:
            StrainMatrix = strain*StrainMode[Modei-1]
            F = crydef.strain2deformgrad(StrainMatrix)
            BaseVecNew = BaseVec.dot(F)
            crylat.print_lattice(BaseVecNew) 
            eglue.write_base_vec(BaseName, BaseVecNew)
            if ParaIn['StableEst']:
                eglue.update_input_for_geomopt(BaseName, FlagEG='E')
            else:
                eglue.update_input_for_geomopt(BaseName, FlagEG='EO')
        elif Modei == 0:
            #Single Point
            eglue.update_input_for_geomopt(BaseName, FlagEG="EO")
        elif Modei == -1:
            #Full geom optimization
            eglue.update_input_for_geomopt(BaseName, FlagEG=Geomi, KeepSym=True)
        time_start = esutils.GetNowTime()
        if ParaIn['EnergyCode'].lower() in RUNDIRECTCODE:
            eglue.run(NP=1, BaseName=BaseName)
        else:
            os.system(RunStr + ">> FP_OUT")        
        time_end = esutils.GetNowTime()
        logger.info('Time cost: {}'.format(time_end - time_start))
        if RemLargFile:
            eglue.rem_files(BaseName)
        E0 = get_org_strain_e(BaseName, flag_se, EnergyCode=ParaIn['EnergyCode'])
        #logger.info(E0)
        np.savetxt(Efilename, E0)
        if Modei == -1:
            eglue.update_structure(BaseName, dstpath="../")
            os.chdir("../")
        elif Modei == 0:
            os.chdir("../../")
    return E0

if __name__ == '__main__':
    logger.info('elastic3rd.elastic is excluded')
else:
    logger.info('elastic3rd.elastic is imported')

################################## END ################################################## 
