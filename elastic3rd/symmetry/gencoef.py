#!python
#
import os
import numpy as np
import elastic3rd.crystal.deform as crydef
import elastic3rd.symmetry.symmetry as essym
import elastic3rd.symmetry.symmetry_stress as esstr

def get_strainmode(CrystalType, Ord, flag_se="s"):
    '''
    Get the build-in strain mode
    Parameter
    ---------
        CrystalType: str
            The crystal type
        Ord: int
            The order of elasti constants
        flag_se: str
            Specify strain-energy (e) or strain-stress (s) method
    Return
    ------
        StrainMode: np.ndarray
            The strain modes for current crystaltype and elastic ord under e/s method
    '''
    coef_path = os.path.dirname(__file__)
    CrystalType = CrystalType.upper()
    pre_filename = coef_path + "/COEFS/" + str(int(Ord)) + "/"
    if flag_se.lower() == "e":
        mid_filename = "ENERGY/"
    elif flag_se.lower() == "s":
        mid_filename = "STRESS/"
    else:
        excpt_str = "No such value(" + flag_se + ") in flag_se"
        raise Exception(excpt_str)
    filedir = pre_filename + mid_filename
    alltype = os.listdir(filedir)
    if Ord == 2:
        if CrystalType in ["C1", "C2", "C"]:
            CrystalType = "C"
        elif CrystalType in ["H1", "H2", "H"]:
            CrystalType = "H"

    if CrystalType in alltype:
        filename = filedir + CrystalType
    else:
        excpt_str = "No such value(" + CrystalType + ") in CrystalType"
        raise Exception(excpt_str)
    StrainMode = read_strainmode(filename)
    return StrainMode

def read_strainmode(filename):
    '''
    Read the build-in strain modes
    Parameter
    ---------
        filename: str
            The filename of the strain mode file
    Return
    ------
        strainmode: np.ndarray
            The strain mode, nx6
    '''
    strainmode_m = np.loadtxt(filename)
    n_strainmode = int(len(strainmode_m)/3)
    strainmode = np.zeros((n_strainmode, 6))
    for i in range(0, n_strainmode):
        strainmode_single = strainmode_m[3*i:3*(i+1)]
        strainmode[i] = crydef.matrix2vec(strainmode_single)
        #print(strainmode)
    return strainmode
