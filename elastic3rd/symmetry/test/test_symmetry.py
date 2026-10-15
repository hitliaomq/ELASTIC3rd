#!python
#
import elastic3rd.symmetry.symmetry as essym
import numpy as np

def test_gen_strain():
    StrainModes = essym.gen_strain(1)
    print(StrainModes)

def test_Num2IJK():
    Num = 216
    N = 6
    Ord = 3
    IJK = essym.Num2IJK(Num, N, Ord)
    print(IJK)

def test_gen_strainord():
    strain = np.array([[1, 0, 0], [0, 1, 1], [0, 1, 0]])
    (StrainOrd, StrainOrdCoef) = essym.gen_strainord(strain)
    print(StrainOrd)
    print(StrainOrdCoef)

def test_gen_strain_coef():
    strain = np.array([[1, 0, 0], [0, 1, 1], [0, 1, 0]])
    (StrainOrd, StrainOrdCoef) = essym.gen_strainord(strain)
    (CoefUniq, CoefUniqCoef) = essym.gen_strain_coef(StrainOrd, StrainOrdCoef, Ord=3)
    print(CoefUniq)
    print(CoefUniqCoef)
    print(type(CoefUniq))

def test_gen_cijk_coef():
    CrystalType = 'c1'
    Strain = np.array([[1, 0, 0], [0, 1, 1], [0, 1, 0]])
    Ord = 3
    CoefStrain = essym.gen_cijk_coef(CrystalType, Strain, Ord)
    print(CoefStrain)
    print(type(CoefStrain))

def test_gen_strain_mode():
    CrystalType = 'c1'
    Ord = 3
    (StrainModeCoef, StrainMode) = essym.gen_strain_mode(CrystalType, Ord)
    print(StrainModeCoef)
    print(StrainMode)
    print(StrainModeCoef.coef2)
    print(StrainModeCoef.coef3)

def test_gen_cijall():
    cijAll = essym.gen_cijall(Ord=2)
    print(cijAll)

def test_cijk2num():
    Cijk = [111, 166]
    CijkNum = essym.cijk2num(Cijk, Ord=3)
    print(CijkNum)

def test_num2cijk():
    num = np.array([1, 8])
    num = 1
    Ord = 3
    cijk = essym.num2cijk(num, Ord)
    print(cijk)

def test_CoefForSingleMode():
    CrystalType = 'c1'
    Strain = np.array([[1, 1, 0, 2, 0, 0]])
    Ord = 3
    (Cijk, StrainModeCoef, StrainMode) = essym.CoefForSingleMode(CrystalType, Ord, Strain)
    print(Cijk.coef2)
    print(Cijk.coef3)
    print(StrainModeCoef)
    print(StrainMode)

def test_print_cijk():
    CrystalType = 'c1'
    Ord = 3
    Cijk = essym.print_cijk(CrystalType, Ord)
    print(Cijk)
    print(type(Cijk))

def test_Convert2indepen():
    Cijk = [44, 66, 13]
    CrystalType = 'c1'
    Ord = 2
    (Cijk1, coef1) = essym.Convert2indepen(Cijk, Ord=Ord, CrystalType=CrystalType)
    assert(Cijk1 == [44, 44, 12])
    assert(coef1 == [1, 1, 1])

    CrystalType = 'h1'
    (Cijk2, coef2) = essym.Convert2indepen(Cijk, Ord=Ord, CrystalType=CrystalType)
    assert(Cijk2 == [44, [11, 12], 13])
    assert(coef2 == [1, [0.5, -0.5], 1])

    Cijk = [222, 166, 113]
    CrystalType = 'c1'
    Ord = 3
    (Cijk, coef3) = essym.Convert2indepen(Cijk, Ord=Ord, CrystalType=CrystalType)
    assert(Cijk == [111, 155, 112])
    assert(coef3 == [1, 1, 1])
    

test_CoefForSingleMode()
