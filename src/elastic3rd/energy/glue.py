#!python
#
import os
import shutil
import numpy as np

def write_energyrun(RunStr):
    EPath = os.path.dirname(__file__)
    FileName = EPath + "\energyrun"
    fopen = open(FileName, 'w')
    fopen.write(RunStr)
    fopen.close()

def run():
    EPath = os.path.dirname(__file__)
    FileName = EPath + "\energyrun"
    fopen = open(FileName, 'r')
    for eachline in fopen:
        eachline = eachline.strip("\n")
        RunStr = eachline.strip()
        print(RunStr)
    fopen.close()

def multi_energy(energy):
    #in unit: eV
    EV2H = 0.0367493
    EV2KJ = 96.4853104
    EV2KC = 23.0605441
    Energy = np.array([energy, EV2H*energy, EV2KJ*energy, EV2KC*energy])
    return Energy


def DeleteFiles(path, remainDirsList, filesList):
    dirsList = []
    dirsList = os.listdir(path)
    for f in dirsList:
        if f not in remainDirsList:
            filePath = os.path.join(path,f)
            if os.path.isdir(filepath):
                shutil.rmtree(filepath, True)
        if f in filesList:
            filepath = os.path.join(path,f)
            os.remove(f)

'''
if __name__ == "__main__":
    path=os.getcwd()+"\\"
    #remainfile
    filesList=['a.txt','b.txt']
    #remain folder
    dirsList=['test']
    DeleteFiles(path,fileList,dirsList)
'''