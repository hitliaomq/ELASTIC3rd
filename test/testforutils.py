#!python
# filename: testforutils.py
# function: test for the esutils.py
#           
# Author: Mingqing Liao
# E-mail: liaomq1900127@163.com
# FGMS @ Harbin Institute of Technology(HIT)

import elastic3rd
import elastic3rd.esutils as esutils

##Test for print_logo
esutils.print_logo()

##Test for read_input
ParaIn = esutils.read_input()
print(ParaIn)

##Test for write_input
FileName = 'INPUT-test'
esutils.write_input(ParaIn, FileName)

##Test for read_strainmode
strainmode = esutils.read_strainmode()
print(strainmode)

##Test for gen_strain_list
StrainList = esutils.gen_strain_list(ParaIn)
print(StrainList)

