#!python
#
import re
import os
import shutil
import linecache
import numpy as np
import elastic3rd.energy.glue as glue

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
    FileName = "OSZICAR"
    fopen = open(FileName, 'r')
    for eachline in fopen:
        linei = eachline.split("=")
        if len(linei) > 2:
            energy = linei[2].strip().split(" ")[0]
    fopen.close()
    energy = float(energy)
    # eV, Hartree, kJ/mol, kcal/mol
    Energy = glue.multi_energy(energy)
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