#!python
# filename: testforsymmetry.py
# function: test for the coefficients(include A2 and A3) generation
#           in symmetry module.
#           
# Author: Mingqing Liao
# E-mail: liaomq1900127@163.com
# FGMS @ Harbin Institute of Technology(HIT)
# 

import numpy as np
import elastic3rd
import elastic3rd.symmetry.symmetry as essym

###########################################################################
#TEST FOR SYMMETRY MODULE
#Two main parts:
#    1. test for symmetry script
#In each test of script, there are two main functions:
#    a. test for a given strain mode
#    b. test for a given symmetry

print("####################Start test for symmetry module#########################")

print("*************Start of the test for symmetry script*************")
print("1. This is the test for a specific strain mode.")
# Generate the coefficients for a given strain mode and specific crystal type
CrystalType = 'c1'
Ord = 3
Strain = np.array([[1, 0, 0, 2, 2, 0]])
(Cijk,StrainModeCoef,StrainMode)=essym.CoefForSingleMode(CrystalType, Ord, Strain)
print(StrainMode)
print("------End of test for a specific strain mode.------\n")

print("2. This is the test for a specific crystal symmetry.")
#Generate the coefficient for specific crystal type
CrystalType = 'c1'
Ord = 3
(StrainModeCoef, StrainMode) = essym.gen_strain_mode(CrystalType, Ord)
print(StrainMode)
print(StrainModeCoef.coef3)
print(StrainModeCoef.coef2)
print("------End of test for a given symmetry.------\n")
print("*************End of test for symmetry script*************")
#End of test for symmetry script

print("####################End test for symmetry module###########################")