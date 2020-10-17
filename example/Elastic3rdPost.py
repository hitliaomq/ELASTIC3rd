#!python
# filename: Elastic3rdPost.py
# function: post for elastic3rd
#           
# Author: Mingqing Liao
# E-mail: liaomq1900127@163.com
# FGMS @ Harbin Institute of Technology(HIT)


import elastic3rd.post.post as espost
import elastic3rd.crystal.lattice as crylat
import elastic3rd.esutils as esutils
from sys import argv
import os

if len(argv) > 1:
    #The first arg must be the INPUT file
    INPUT = argv[1]
else:
    INPUT = "INPUT"
if not os.path.isfile(INPUT):
    raise Exception("ERROR:" + INPUT + "file not exists")

INPUT = "INPUT"
ParaIn = esutils.read_input(INPUT)
eglue = __import__("energy." + ParaIn['EnergyCode'], fromlist = ParaIn['EnergyCode'])

BaseName = ParaIn['BaseName']

BaseVec = eglue.get_base_vec(BaseName)
V0 = crylat.print_lattice(BaseVec)

V0 = 163.193194972   #The volume of undeformed structure
Flag_Fig = 1
Flag_Ord = 3
EEnergy = "EEnergy.txt"

#Post for a given symmetry, no specific the strain mode
(C2, C3) = espost.post(V0, Flag_Fig, Flag_Ord, EEnergy, INPUT)
print(C2)
print(C3)

'''
#Post for some specific the strain mode
STRAINMODE = "STRAINMODE"
(C2, C3) = espost.post_mode(V0, Flag_Fig, Flag_Ord, EEnergy, INPUT, STRAINMODE)
print(C2)
print(C3)
'''