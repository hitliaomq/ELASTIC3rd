#!python
#
import re
import os
import shutil
import linecache
import numpy as np

def get_base_vec(BaseName):
    BaseVec = np.zeros((3, 3))
    FileName = "POSCAR"
    fopen = open(FileName, 'r')
    count = 0
    for eachline in fopen:
        if count > 1:
            linei = eachline.strip('\n').strip()
            linei = re.split('\s+', linei)
            for j in range(0, 3):
                BaseVec[count - 2][j] = float(linei[j])
        if count > 3:
            return BaseVec
        count = count + 1
    return BaseVec

def write_base_vec(BaseName, BaseVec):
    FileName = "POSCAR"
    fopen = open(FileName, 'r')
    tmpopen = open('tmpfile', 'a')
    count = 1
    lines = []
    for eachline in fopen:
        if (count > 2) and (count < 6):
            for j in range(0, 3):
                tmpopen.write("      ")
                tmpopen.write("%.15f" % BaseVec[count - 3][j])
            tmpopen.write("\n")
        else:
            tmpopen.write(eachline)
        count = count + 1
    fopen.close()
    tmpopen.close()
    os.remove(FileName)
    os.rename('tmpfile', FileName)

def run(NP, BaseName):
    RunStr = "mpirun -np " + str(int(NP)) + " vasp_std"
    return RunStr

def get_energy(BaseName):
    EV2H = 0.0367493
    EV2KJ = 96.4853104
    EV2KC = 23.0605441
    FileName = "OSZICAR"
    fopen = open(FileName, 'r')
    for eachline in fopen:
        linei = eachline.split("=")
        if len(linei) > 2:
            energy = linei[2].strip().split(" ")[0]
    fopen.close()
    energy = float(energy)
    # eV, Hartree, kJ/mol, kcal/mol
    Energy = np.array([energy, EV2H*energy, EV2KJ*energy, EV2KC*energy])
    return Energy

def copy_files(BaseName, Path):
    #TODO: need to be improved
    shutil.copyfile("INCAR", Path + "/" + "INCAR")
    shutil.copyfile("POSCAR", Path + "/" + "POSCAR")
    shutil.copyfile("POTCAR", Path + "/" + "POTCAR")
    shutil.copyfile("KPOINTS", Path + "/" + "KPOINTS")

'''
#TEST
os.mkdir("Al")
copy_files("Al", "Al")
'''