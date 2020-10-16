#!python
#
import numpy as np
import re
import os

def read_input(INPUT = "INPUT"):
    '''
    Read the input file and set default value for Input parameters
    Parameters
    ----------
        INPUT: str
            The file name of the input file
                The key and value in the INPUT file should be divided by space
                The key is capital sensitive and the value is not
    Return
    ------
        ParaIn: dict
            The dict of input parameters
    '''
    ParaIn = {'Continue': 0, 'EnergyCode' : 'CASTEP', 'EnergyRun' : 1, 'CrystalType': 'Cubic1', 'Ord': 3, 
        'MaxStrain' : 5, 'NP' : 1, 'STEPS' : 5, 'BaseName' : 'Si', 'CalMode': 'W', 'FlagSE': 'E'}
    fopen = open(INPUT, 'r')
    for eachline in fopen:
        eachline = eachline.strip('\n')
        linei = re.split('\s+', eachline)
        if len(linei) >= 2 :
            flag = linei[0]
            value = linei[1]
            try:
                ParaIn[flag] = float(value)
            except:
                ParaIn[flag] = value
            #ParaIn[flag] = value
    ParaIn['EnergyCode'] = ParaIn['EnergyCode'].lower()
    if ParaIn['FlagSE'].lower() == 's':
        print("Warning: The FlagSE keyword is under developing, and not work. It will set FlagSE to E.")
        ParaIn['FlagSE'] = 'E'
    fopen.close()
    return ParaIn

def write_input(ParaIn, FileName = "INPUT"):
    '''
    Write INPUT file according ParaIn dict
    Parameters
    ----------
        ParaIn: dict
            The dict of INPUT parameters
        FileName: str
            The file name of the INPUT file
    Return
    ------
        None, but it will write out an input file
    '''
    with open(FileName, 'w+') as fopen:
        for key in ParaIn:
            fopen.write('{}  {}\n'.format(key, ParaIn[key]))

def read_strainmode(STRAINMODE = "STRAINMODE"):
    '''
    Read strain mode file
    Parameters
    ----------
        STRAINMODE: str
            The filename of strain mode file, STRAINMODE by default
    Return
    ------
        strainmode: np.ndarray
            strain modes, 2D matrix, e.g. [[1, 2, -0.5, 1, 1, 0]]
    '''
    strainmodetmp = np.loadtxt(STRAINMODE)
    m = strainmodetmp.shape
    if len(m) == 1:
        strainmode = np.zeros((1, 6))
        strainmode[0, :] = strainmodetmp
    else:
        strainmode = strainmodetmp
    return strainmode

def print_parain(ParaIn):
    '''
    Print the input parameters
    Parameters
    ----------
        ParaIn: dict
            The dict of input parameters
    Return
    ------
        None
    '''
    if ParaIn['Continue']:
        RunMode = "CONTINUE"
    else:
        RunMode = "NEW"
    if ParaIn['FlagSE'].lower() == 'e':
        Method = "Strain-Energy method"
    else:
        Method = "Strain-Stress method"    
    print("The Runing mode is : " + RunMode)
    print("The BaseName is : " + ParaIn['BaseName'])
    print("The Method used here is: " + Method)
    if ParaIn['CalMode'].lower() == 's':
        print("The strain mode(s) is(are) specified.")
    print("The crystal symmetry is : " + ParaIn['CrystalType'])
    print("The order of elastic is : " + str(ParaIn['Ord']))
    print("The code for calculating the energy : " + ParaIn['EnergyCode'].upper())
    print("The maximum of the strain : " + str(ParaIn['MaxStrain']))
    print("The steps run in every mode : " + str(ParaIn['STEPS']))
    print("The core used in energy calculate : " + str(ParaIn['NP']) + "\n")
    #print("\n")

def creat_folder(path):
    '''
    Create folder if not exist, if exist, skip
    Parameters
    ----------
        path: str
            The path of the folder
    Return
    ------
        None
    '''
    if os.path.exists(path):
        print("Warning: The folder(" + path + ") exists!\n")
    else:
        os.mkdir(path)

def print_e(E):
    '''
    Print the energy
    Parameters
    ----------
        E: list
            The energy in different unit, it is generated by elastic3rd.energy.glue.multi_energy function
    Return
    ------
        None
    '''
    print("       eV             Hartree            kJ/mol            kcal/mol")
    E_str = []
    for i in range(0, 4):
        strparam = "%.10f"
        E_str.append(strparam % E[i])
    print("  ".join(E_str))

def gen_strain_list(ParaIn):
    '''
    Generate the strain list according to the input parameters
        The range of the strain is -MaxStrain to MaxStrain
        The number of point is determined by the STEP
    Parameters
    ----------
        ParaIn: dict
            The dict of the input parameters
    Return
    ------
        StrainList: np.ndarray
            The list of the strain
    '''
    Step = int(ParaIn['STEPS'])
    MaxStrain = float(ParaIn['MaxStrain'])
    if Step % 2 == 0:
        Step = Step + 1
    if Step < 3:
        print("ERROR:The STEPS is too small, at least 3!")
        return
    StepStrain = 2*MaxStrain/(Step - 1)
    EndStrain = MaxStrain + StepStrain
    StrainList = np.arange(-MaxStrain, EndStrain, StepStrain)
    return StrainList

def iscontinue(Path, ModeOrStrain, flag_se = "e"):
    '''
    If current mode or strain is calculated or not
        It is judged by the existance of the result file
            In the mode, if current mode is finished, there should exists a file named Energy_Mode.txt
            In the strain, if current strain is finished, there should exists a file named Energy_Strain.txt
    Parameters
    ----------
        Path: str
            The path of the mode or strain
        ModeOrStrain: str
            Mode or Strain
        flag_se: str
            e or s. e for strain-energy method and s for strain-stress method.
                Currently, the 's' is not work
    Return
    ------
        flag: 0 or 1
            If current mode or strain is calculated, return 1, orhterwise return 0
    '''
    if flag_se == "e":
        filename = Path + "/" + "Energy_" + ModeOrStrain + ".txt"
    elif flag_se == "s":
        filename = Path + "/" + "Stress_" + ModeOrStrain + ".txt"
    if os.path.exists(filename):
        flag = 1
    else:
        flag = 0
    return flag

def print_logo():
    '''
    Pring the LOGO
    '''
    print("+=============================================================+")
    print("| EEEEE L       A   SSSSS TTTTT IIIII CCCCC 33333 RRRRR DDD   |")
    print("| E     L      A A  S       T     I   C         3 R   R D  D  |")
    print("| EEEE  L     A   A SSSSS   T     I   C      3333 RRRRR D   D |")
    print("| E     L     AAAAA     S   T     I   C         3 R  R  D  D  |")
    print("| EEEEE LLLLL A   A SSSSS   T   IIIII CCCCC 33333 R   R DDD   |")
    print("+-------------------------------------------------------------+")
    print("|             Version: ELASTIC3RD  version 2.5.0              |")
    print("|                Date: 2019-05-26                             |")
    print("|              Author: Liao Mingqing                          |")
    print("|              E_mail: liaomq1900127@163.com                  |")
    print("|     FGMS @ Harbin Institute of Technology(HIT)              |")
    print("|    URL:https://github.com/hitliaomq/ELASTIC3rd              |")
    print("|Please Cite:                                                 |")
    print("| Mingqing Liao,ELASTIC3rd,(2018).doi:10.5281/zenodo.2525580  |")
    print("+=============================================================+")


