#!python
import re
import os
import numpy as np
from sys import argv
import shutil

def get_pos(BaseName):
    #Stress = np.zeros((3, 3))
    #
    castepfile = BaseName + ".castep"
    fopen = open(castepfile, 'r')
    position_tmp = []
    flag = 0
    #flag_be  the flag of begin or end of the position
    flag_be = 0  
    count = 0
    i = 0
    for eachline in fopen:
        eachline = eachline.strip('\n')
        eachline = eachline.strip();
        #eachlinei = eachline.split("*")
        if 'Cell Contents' in eachline:
            flag = 1
        if flag :
            if 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' in eachline:
                if flag_be == 0:
                    flag_be = 1
                else:
                    flag_be = 0
                    count = 0
                    flag = 0
                    position_tmp = []
            if flag_be == 1:
                count = count + 1
        if count < 5:
            continue
        else:
            pos = re.split(r'\s+', eachline)
            #print(pos)
            posi = [pos[1]]
            posi.extend(pos[3:6])
            position_tmp.append(posi)
            position = position_tmp
            #Stress[i, :] = Stressi[3:5]
            
            i = i + 1
            #if count > 9:
            #    break
    fopen.close()
    #print(position)
    #if flag == 0:
    #    print "Please set calculate_stress as ture."
    return position

def update_pos(BaseName, NewPos):
    FileName = BaseName + ".cell"
    fopen = open(FileName, 'r')
    tmpopen = open('tmpfile', 'a')
    flag = 0
    count = 0
    lines = []
    for eachline in fopen:
        linei = eachline.strip('\n')        
        linei = linei.strip();
        #print eachline
        if linei == "%BLOCK POSITIONS_FRAC":
            tmpopen.write("%BLOCK POSITIONS_FRAC\n")
            flag = 1
            continue
        elif linei == "%ENDBLOCK POSITIONS_FRAC":
            flag = 0
        if flag:
            for j in range(0,4):
                tmpopen.write("    ")
                tmpopen.write("%s" % NewPos[count][j])
            tmpopen.write("\n")
            count = count + 1
        else:
            tmpopen.write(eachline)
    fopen.close()
    tmpopen.close()
    os.rename(FileName, BaseName + "-org.cell")
    #os.remove(FileName)
    os.rename('tmpfile', FileName)

def get_basevec(BaseName):
    BaseVec = np.zeros((3, 3))
    FileName = BaseName + ".castep"
    fopen = open(FileName, 'r')
    flag = 0
    count = 1
    i = 0
    for eachline in fopen:
        eachline = eachline.strip('\n')
        eachline = eachline.strip();
        #eachlinei = eachline.split("*")
        if 'Unit Cell' in eachline:
            flag = 1
        if flag :
            count = count + 1
            if count < 5:
                continue
            else:
                BaseVeci = re.split(r'\s+', eachline)
                BaseVec[i, :] = BaseVeci[0:3]
                i = i + 1
                if count > 6:
                    flag = 0
                    count = 1
                    i = 0
    fopen.close()
    return BaseVec

def update_basevec(BaseName, BaseVec):
    FileName = BaseName + ".cell"
    fopen = open(FileName, 'r')
    tmpopen = open('tmpfile', 'a')
    flag = 0
    count = 0
    lines = []
    for eachline in fopen:
        linei = eachline.strip('\n')        
        linei = linei.strip();
        #print eachline
        if linei == "%BLOCK LATTICE_CART":
            tmpopen.write("%BLOCK LATTICE_CART\n")
            flag = 1
            continue
        elif linei == "%ENDBLOCK LATTICE_CART":
            flag = 0
        if flag:
            for j in range(0,3):
                tmpopen.write("      ")
                tmpopen.write("%.15f" % BaseVec[count][j])
            tmpopen.write("\n")
            count = count + 1
        else:
            tmpopen.write(eachline)
    fopen.close()
    tmpopen.close()
    os.remove(FileName)
    os.rename('tmpfile', FileName)

#BaseName = argv[1]
BaseName = 'C'   #for test
shutil.copyfile(BaseName + ".cell", BaseName + "-bk.cell")
position = get_pos(BaseName)
print(position)
update_pos(BaseName, position)
#BaseVec = get_basevec(BaseName)
#update_basevec(BaseName, BaseVec)