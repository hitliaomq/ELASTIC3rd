#!python
#
import numpy as np
import re
import os

def read_input(INPUT = "INPUT"):
    #Read the input
    #the default value
    ParaIn = {'EnergyCode' : 'CASTEP', 'EnergyRun' : 1, 'MaxStrain' : 5, 'NP' : 1, 'STEPS' : 5, 'BaseName' : 'Si'}
    fopen = open(INPUT, 'r')
    for eachline in fopen:
        eachline = eachline.strip('\n')
        linei = re.split('\s+', eachline)
        if len(linei) >= 2 :
            flag = linei[0]
            value = linei[1]
            ParaIn[flag] = value
    ParaIn['EnergyCode'] = ParaIn['EnergyCode'].lower()
    fopen.close()
    return ParaIn

def print_parain(ParaIn):
    print("The BaseName is : " + ParaIn['BaseName'])
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

def print_logo():
    print("+==========================================================+")
    print("|  EEEEE L     SSSSS TTTTT IIIII CCCCC 33333 PPPPP Y   Y   |")
    print("|  E     L     S       T     I   C         3 P   P  Y Y    |")
    print("|  EEEE  L     SSSSS   T     I   C      3333 PPPPP   Y     |")
    print("|  E     L         S   T     I   C         3 P       Y     |")
    print("|  EEEEE LLLLL SSSSS   T   IIIII CCCCC 33333 P       Y     |")
    print("+----------------------------------------------------------+")
    print("|           Version: ELASTIC3PY                            |")
    print("|              Date: 2018-12-22                            |")
    print("|            Author: Liao Mingqing                         |")
    print("|            E_mail: liaomq1900127@163.com                 |")
    print("|   FGMS @ Harbin Institute of Technology(HIT)             |")
    print("|  URL:https://github.com/hitliaomq/ELASTIC3rd             |")
    print("|Please Cite:                                              |")
    print("|Mingqing Liao,ELASTIC3rd,(2018).doi:10.5281/zenodo.2525581|")
    print("+==========================================================+")

'''
#print_logo()
ParaIn = read_input('INPUT')
print ParaIn
StrainList = gen_strain_list(ParaIn)
print StrainList
'''