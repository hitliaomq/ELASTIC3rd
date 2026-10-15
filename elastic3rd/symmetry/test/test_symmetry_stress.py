#!python

import numpy as np
import elastic3rd.symmetry.symmetry_stress as symstr
import elastic3rd.symmetry.symmetry as essym

def test_cauchy2pk():
    LagStrain = 0.05*np.array([[1, 0, 0], [0, 1, 1], [0, 1, 0]])
    Cauchy = np.array([[15, 0, 0], [0, 15, 0], [0, 0, 0]])
    (PK, PKvec) = symstr.cauchy2pk(Cauchy, LagStrain)
    print(PK)
    print(PKvec)

def test_gen_strain_coef():
    strain = np.array([[1, 0, 0], [0, 1, 1], [0, 1, 1]])
    (StrainOrd, StrainOrdCoef) = essym.gen_strainord(strain)
    (CoefUniq, CoefUniqCoef) = symstr.gen_strain_coef(StrainOrd, StrainOrdCoef, Ord=3)
    print(StrainOrd)
    print(StrainOrdCoef)
    print(CoefUniq)
    print(CoefUniqCoef)
    print(type(CoefUniq))

def test_gen_cijk_coef():
    CrystalType = 'c1'
    Strain = np.array([[1, 0, 0], [0, 1, 1], [0, 1, 0]])
    Ord = 3
    CoefStrain = symstr.gen_cijk_coef(CrystalType, Strain, Ord)
    print(CoefStrain)

def test_CoefForSingleMode():
    CrystalType = 'c1'
    Strain = np.array([[1, 1, 1, 0, 0, 0]])
    Strain = np.array([[1, 1, 0, 2, 4, 0]])
    Ord = 3
    (Cijk, StrainModeCoef, StrainMode) = symstr.CoefForSingleMode(CrystalType, Ord, Strain)
    #print(Cijk.coef2)
    print(Cijk.coef3)
    print(StrainModeCoef.coef2)
    print(StrainMode)

def test_gen_fullrank_db():
    CrystalType = 'c1'
    Ord = 3
    StrainMode = symstr.gen_fullrank_db(CrystalType, Ord)
    print(StrainMode)

def test_gen_strain_mode():
    CrystalType = 'c1'
    Ord = 3
    (StrainModeCoef, StrainMode) = symstr.gen_strain_mode(CrystalType, Ord)
    print(StrainModeCoef)
    print(StrainMode)
    print(StrainModeCoef.coef2)
    print(StrainModeCoef.coef3)

#test_gen_strain_mode()

test_CoefForSingleMode()