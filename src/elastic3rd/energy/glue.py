#!python
#
import os
import shutil
import numpy as np

def write_energyrun(RunStr):
    '''
    Write string to energyrun file
        This function is used to write correct command for calling first principles code
    Parameter
    ---------
        RunStr: str
            The string for run first principles code, e.g. 'mpirun -n 8 vasp_std'
    Return
    ------
        None
    '''
    EPath = os.path.dirname(__file__)
    FileName = os.path.join(EPath + "energyrun")
    fopen = open(FileName, 'w')
    fopen.write(RunStr)
    fopen.close()

def run():
    '''
    Get the run string from energyrun file
    Parameter
    ---------
        None
    Return
    ------
        RunStr: str
            The command calling first principles code.
    '''
    EPath = os.path.dirname(__file__)
    FileName = os.path.join(EPath + "energyrun")
    fopen = open(FileName, 'r')
    for eachline in fopen:
        eachline = eachline.strip("\n")
        RunStr = eachline.strip()
        print(RunStr)
    fopen.close()
    return RunStr

def multi_energy(energy):
    '''
    Transform the energy from ev to different unit
    Parameter
    ---------
        energy: float
            energy in eV
    Return
    ------
        Energy: np.ndarray
            energy in different unit, in order of Hatree, KJ/mol and KCal/mol
    '''
    EV2H = 0.0367493
    EV2KJ = 96.4853104
    EV2KC = 23.0605441
    Energy = np.array([energy, EV2H*energy, EV2KJ*energy, EV2KC*energy])
    return Energy


def DeleteFiles(path, remainDirsList, filesList):
    '''
    Detele files (Not used, may be os.walk is a better way)
    Parameter
    ---------
        path: path-like string
            The path to delete files
        remainDirsList: list
        filesList
    '''
    dirsList = []
    dirsList = os.listdir(path)
    for f in dirsList:
        if f not in remainDirsList:
            filepath = os.path.join(path,f)
            if os.path.isdir(filepath):
                shutil.rmtree(filepath, True)
        if f in filesList:
            filepath = os.path.join(path,f)
            os.remove(f)

