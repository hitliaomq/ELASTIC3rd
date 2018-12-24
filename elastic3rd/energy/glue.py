#!python
#
import os

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
        print RunStr
    fopen.close()
