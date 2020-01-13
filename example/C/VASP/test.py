#!/home/msi/.pyenv/versions/2.7.16/bin/python2.7
#
import os
import sys
import numpy as np

CurPath = os.path.dirname(__file__)
CurPath = CurPath.split('\\')
#print CurPath
ParPath = "\\".join(CurPath[0:-2])
if not ParPath in sys.path:
    sys.path.append(ParPath)

import elastic3rd
import elastic3rd.elastic
import elastic3rd.esutils

elastic3rd.elastic.elastic3()

#import elastic3rd.esutils

#elastic3rd.esutils.print_logo()

'''
ParaIn = elastic3rd.esutils.read_input("INPUT")
BaseName = ParaIn['BaseName']
dstpath = BaseName + "/Mode0"
print dstpath
if ParaIn['Continue']:
    flag_continue = elastic3rd.esutils.iscontinue(dstpath, "Mode")
    print flag_continue

filename = dstpath + "/" + "Energy_Mode.txt"
if os.path.exists(filename):
    flag = 1
    print flag
'''