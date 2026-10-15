#!python
#
import numpy as np
import os
import re

def sortByOrder(BaseName, flag, OrdMax = 7, OrdMin = 3):
    resultfolder = BaseName + "-" + flag + "-Result"
    PreName = resultfolder + "/" + BaseName + "-" + flag + "-"
    Ext = ".txt"
    for k in ('C2', 'C3'):
        C0 = np.loadtxt(PreName + k + "-Ord" + str(OrdMin) + Ext)
        Nrow = len(C0)
        if k == 'C2':
            N = 4
        elif k == 'C3':
            N = 7
        else:
            raise Exception
        for i in range(1, N):
            Csave = np.zeros([Nrow, OrdMax + 2 - OrdMin])
            Csave[:, 0] = C0[:, 0]
            for j in range(OrdMin, OrdMax + 1):
                filename = PreName + k + "-Ord" + str(j) + Ext
                C = np.loadtxt(filename)
                Csave[:, 1 + j - OrdMin] = C[:, i]
            #Write out as BaseName-flag-C2(3)-NO.
            np.savetxt(PreName + k + "-" + str(i) + Ext, Csave)

def check_input(input_sets = {'MaxStrain' : 5.00, 'MinStrain' : -5.00}):
    #Check the INPUT file
    FileName = 'INPUT.e3rd'
    #Params = {'task': 'GeometryOptimization'}
    fopen = open(FileName, 'r')
    tmpopen = open('tmpfile', 'a')
    for eachline in fopen:
        eachlinei = eachline.strip('\n')
        #linei = eachlinei.split(" ")
        linei = re.split(r'\s+', eachline)
        if len(linei) >= 2 :
            flag = linei[0]
            value = linei[1]
            if flag in input_sets:
                if value == input_sets[flag]:
                    tmpopen.write(eachline)
                else:
                    strout = flag + "    " + str(input_sets[flag]) + "\n" 
                    tmpopen.write(strout)
                input_sets.pop(flag)
            else:
                tmpopen.write(eachline)
    fopen.close()
    if input_sets:
        for key in input_sets.keys():
            strout = key + "    " + str(input_sets[key]) + "\n"
            #print(strout)
            tmpopen.write(strout)
    tmpopen.close()
    os.remove(FileName)
    os.rename('tmpfile', FileName)
