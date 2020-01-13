#!python
#
import numpy as np
import re
import os

def read_input(INPUT = "INPUT"):
    #Read the input
    #the default value
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
    #FileName = "INPUT"
    fopen = open(FileName, 'w')
    fopen.write("BaseName        %s\n" % ParaIn['BaseName'] )
    fopen.write("CalMode         %s\n" % ParaIn['CalMode'])
    fopen.write("FlagSE          %s\n" % ParaIn['FlagSE'])
    fopen.write("Continue        %d\n" % ParaIn['Continue'])
    fopen.write("CrystalType     %s\n" % ParaIn['CrystalType'])
    fopen.write("Ord             %d\n" % ParaIn['Ord'])
    fopen.write("EnergyCode      %s\n" % ParaIn['EnergyCode'])
    fopen.write("EnergyRun       %d\n" % ParaIn['EnergyRun'])
    fopen.write("MaxStrain       %.2f\n" % ParaIn['MaxStrain'])
    fopen.write("STEPS           %d\n" % ParaIn['STEPS'])
    fopen.write("NP              %d\n" % ParaIn['NP'])
    fopen.close()

def read_strainmode(STRAINMODE = "STRAINMODE"):
    strainmodetmp = np.loadtxt(STRAINMODE)
    m = strainmodetmp.shape
    if len(m) == 1:
        strainmode = np.zeros((1, 6))
        strainmode[0, :] = strainmodetmp
    else:
        strainmode = strainmodetmp
    return strainmode

def print_parain(ParaIn):
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
    if os.path.exists(path):
        print("Warning: The folder(" + path + ") exists!\n")
    else:
        os.mkdir(path)

def print_e(E):
    print("       eV             Hartree            kJ/mol            kcal/mol")
    E_str = []
    for i in range(0, 4):
        strparam = "%.10f"
        E_str.append(strparam % E[i])
    print("  ".join(E_str))

def gen_strain_list(ParaIn):
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
    #In mode step, the energy file named as Energy_Mode
    #In strain step, the energy file named as Energy_Strain
    #ModeOrStrain equal to Mode or Strain
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
    print("+=============================================================+")
    print("| EEEEE L       A   SSSSS TTTTT IIIII CCCCC 33333 RRRRR DDD   |")
    print("| E     L      A A  S       T     I   C         3 R   R D  D  |")
    print("| EEEE  L     A   A SSSSS   T     I   C      3333 RRRRR D   D |")
    print("| E     L     AAAAA     S   T     I   C         3 R  R  D  D  |")
    print("| EEEEE LLLLL A   A SSSSS   T   IIIII CCCCC 33333 R   R DDD   |")
    print("+-------------------------------------------------------------+")
    print("|             Version: ELASTIC3PY  version 2.4.2              |")
    print("|                Date: 2019-05-26                             |")
    print("|              Author: Liao Mingqing                          |")
    print("|              E_mail: liaomq1900127@163.com                  |")
    print("|     FGMS @ Harbin Institute of Technology(HIT)              |")
    print("|    URL:https://github.com/hitliaomq/ELASTIC3rd              |")
    print("|Please Cite:                                                 |")
    print("| Mingqing Liao,ELASTIC3rd,(2018).doi:10.5281/zenodo.2525580  |")
    print("+=============================================================+")

'''
#print_logo()
ParaIn = read_input('INPUT')
print ParaIn
StrainList = gen_strain_list(ParaIn)
print StrainList
'''
