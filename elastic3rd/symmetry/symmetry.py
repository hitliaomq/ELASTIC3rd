#!python
#

import numpy as np

def gen_strain():
    pass

def gen_coef(strain):
    pass

def voigtmap(i, j):
    if i * j == 1:
        voigt = 1
    elif i * j == 4:
        voigt = 2
    elif i * j == 9:
        voigt = 3
    elif i * j == 6:
        voigt = 4
    elif i * j == 3:
        voigt = 5
    elif i * j == 2:
        voigt = 6
    else:
        voigt = 0
        print("ERROR: Both i and j shold range from 1 to 3.")
    return voigt
