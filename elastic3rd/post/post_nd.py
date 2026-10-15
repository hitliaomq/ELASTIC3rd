#!python
#
#Post module for numerical differetiation method for general SEM and SSM
#Features:

import elastic3rd.tools.nd_coef as nd_coef
import elastic3rd.post.post2 as espost
import elastic3rd.symmetry.symmetry as essym
import elastic3rd.symmetry.symmetry_stress as symstr
import elastic3rd.esutils as esutils
import warnings
import numpy as np

def numerr_round(a, n=1, k=8):
    '''
    This function is used to round the float caused by numerical error.
        This function will neglect continue numbers(contraled by k) 0 or 9
        e.g numerr_round(2.499999996) = 2.499999996 and numerr_round(2.4999999996) = 2.5
            numerr_round(2.500000004) = 2.500000004 and numerr_round(2.5000000004) = 2.5
            a = 0.135/3
            print(a)
            print(numerr_round(a))
    '''
    if abs(a) < pow(10, -k):
        return 0
    b1 = round(a, n)
    b2 = round(a, n+k)
    if b1 == 0.:
        return numerr_round(a, n=n+1)
    else:
        if b1 == b2:
            return b1
        else:
            return numerr_round(a, n=n+1)

def check_nd_strain(s, e, acc=2):
    '''
    This function is used to check the strain if it is suitable for numerical differentiation method,
    it should follow the following rules:
        1. The strain should be sorted, and corresponding energy/stress
        2. The number of strains should equal acc
        3. The strain should be symmetrical around 0 (optional, if voliet, give an warning)

    Test
    ----
        s = [-2, -3, -1., -4, 1e-10, 1, 2, 3, 4]
        e = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        (s, e) = check_nd_strain(s, e, acc=7)
        print(s)
        print(e)
    '''
    s = [numerr_round(si) for si in s]
    e = sort_x_by_y(e, s)
    s = sorted(s)
    s_no0 = [si for si in s if abs(si) > 1e-8]
    if len(s_no0) < acc:
        raise ValueError('The number of strain (excluding 0) should no less than acc')
    for i in range(int(len(s)/2)):
        if abs(s[i]) != abs(s[-(i+1)]):
            warnings.warn('The strain is not symmetrical around 0')
    return (s, e)

def check_strain0(s, e, P0=0.):
    '''
    Get the energy or strain of inital structure
    '''
    index0 = [i for i, si in enumerate(s) if abs(si) < 1e-8]
    if len(index0) == 0:
        e0 = P0
    else:
        e0 = e.pop(index0[0])
        s.pop(index0[0])
    e = sort_x_by_y(e, s)
    s = sorted(s)
    return (s, e, e0)

def sort_x_by_y(x, y):
    """
    Sort a list of x in the order of sorting y
    **This function is taken from dfttk**
    """
    return [xx for _, xx in sorted(zip(y, x), key=lambda pair: pair[0])]

def get_strain_by_acc(s, e, acc=2, P0=0):
    '''
    Get the strain list according to accuracy, it will take the strain by the following rules:
        1. The eq structure is always in the list (If it exists in s, using that, if not, using P0)
        2. The others are taken from the first acc/2 and the last acc/2. e.g. for acc=2, it will take the first one and the last one
    '''
    (s_no0, e_no0, e0) = check_strain0(s, e, P0=P0)
    half_acc = int(acc/2)
    s_sel = [0] * (half_acc*2 + 1)
    e_sel = [e0] * (half_acc*2 + 1)
    for i in range(half_acc):
        s_sel[i] = s_no0[i]
        s_sel[-(i+1)] = s_no0[-(i+1)]
        e_sel[i] = e_no0[i]
        e_sel[-(i+1)] = e_no0[-(i+1)]
    return (s_sel, e_sel)

#eVpmol2GPa = 160.21719175
def get_coef_nd(s, e, flag_se='s', P0=0, acc=2, EOrd=3):
    '''
    This function is used to get the coefficients (A2 and A3) by numerical differentiation method
    Parameter
    ---------
        s: list/np.ndarray
            The strain list. Note the large strain must be the int times of the minimum strain(absolute value)
        e: list/np.ndarray
            energy density or stress
        flag_se: str
            strain-stress method (s) or strain-energy method (e)
        P0: float/list
            The stress or energy at eq structure
        acc: int
            The accuracy, 2 or 4 or 6 or 8
        EOrd: int (unused now, only for 3 order)
            The order of elastic constant
    Return
    ------
        coef: np.ndarray
            A2 and A3
    '''
    flag_se = flag_se.lower()
    
    e = np.array(e)
    (m, n) = e.shape
    coef_fit = np.zeros((n, 2))

    for i in range(0, n):
        ei = e[:, i]
        si = [sj for sj in s]
        #e0 = e[n_d][i]  
        (si, ei) = check_nd_strain(si, ei, acc=acc)
        (si, ei) = get_strain_by_acc(si, ei, acc=acc, P0=P0)
        
        si = np.array(si)
        ei = np.array(ei)

        if flag_se == "e":
            ndcoef2 = nd_coef.gen_nd_coef_e(si, deg=2)
            ndcoef3 = nd_coef.gen_nd_coef_e(si, deg=3)
        elif flag_se == "s":
            ndcoef2 = nd_coef.gen_nd_coef(si, deg=1)
            ndcoef3 = nd_coef.gen_nd_coef(si, deg=2)
        else:
            raise ValueError('Currently only support flag_se(e and s)')

        coef_fit[i, 0] = np.dot(ei, ndcoef2)
        coef_fit[i, 1] = np.dot(ei, ndcoef3)

    return coef_fit

def post_single(x, E, StrainIn, INPUT="INPUT.e3rd", P0=0, acc=2, EOrd=3):
    (CrystalType, Ord, flag_se, StrainList) = espost.get_post_param(INPUT)
    if flag_se == "e":
        Cij_mode, coef_e, StrainMode = essym.CoefForSingleMode(CrystalType, Ord, StrainIn)
    elif flag_se == "s":
        Cij_mode, coef_e, StrainMode = symstr.CoefForSingleMode(CrystalType, Ord, StrainIn)
    else:
        raise ValueError

    coef_fit = get_coef_nd(x, E, flag_se=flag_se, P0=P0, acc=acc, EOrd=EOrd)
    coef2 = coef_e.coef2
    coef3 = coef_e.coef3
    (C2, C3) = espost.get_cij(coef_fit, coef2, coef3, flag_se=flag_se)

    return (C2, C3)


