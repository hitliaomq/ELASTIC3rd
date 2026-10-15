#!python
#
#symopt: symmetry operater
#
#NOTE: This script is not used now. It's for n-th order elastic constants usage
#
import elastic3rd.symmetry.symmetry as essym
from itertools import combinations, product
import numpy as np
import copy
from scipy.linalg import lu
import math

def symoptgendb_2d(symkey):
    #symoptdb_var = {'2' : [3], '2mm' : [6, 3], '4mm' : [8, 6], '6mm' : [12, 6]}
    symoptdb_var = {'1' : [1], '2' : [3], 'm': [6], '2mm' : [6, 3], 
                '4' : [8], '4mm' : [8, 6], '3' : [10], '3m' : [10, 6], 
                '6' : [12], '6mm' : [12, 6]}
    optnum = symoptdb_var[symkey]
    return optnum

def symoptdb_2d(NO):
    #Note : this is only the basic transformation matrix, it need to be expand
    sqrt3 = math.sqrt(3.)/2.0
    '''
    symdb = {
        3 : np.array([[-1, 0], [0, -1]]), 
        6 : np.array([[-1, 0], [0, 1]]),
        8 : np.array([[0, 1], [-1, 0]]),
        12 : np.array([[0.5, sqrt3], [-sqrt3, 0.5]])
    }
    '''
    symdb = {
        1 : np.array([[1, 0], [0, 1]]), 
        3 : np.array([[-1, 0], [0, -1]]), 
        6 : np.array([[-1, 0], [0, 1]]),
        8 : np.array([[0, 1], [-1, 0]]),
        10 : np.array([[-0.5, sqrt3], [-sqrt3, -0.5]]),
        12 : np.array([[0.5, sqrt3], [-sqrt3, 0.5]]), 
    }
    optmat = symdb[NO]
    return optmat

def cij_voigt(IJK):
    '''
    Convert the tensor index to voigt's notation
        [[1, 1, 1, 2]] -> [16]
    '''
    voigts = []
    try:
        n = len(IJK[0])
    except Exception as e:
        n = len(IJK)
        IJK = [IJK]
    for cij_ten in IJK:
        voigti = []
        for i in range(0, int(n/2)):
            voigt = essym.voigtmap(cij_ten[2*i], cij_ten[2*i+1])
            voigti.append(voigt)
        voigti = sorted(voigti)
        voigti_sum = 0
        for i in range(0, int(n/2)):
            voigti_sum = voigti_sum + voigti[i] * (10 ** (int(n/2) - 1 - i))
        voigts.append(voigti_sum)
    return voigts

def get_cij_uniq_from_ord(Ord, N=3):
    N_total = N ** Ord
    IJK_list = []
    for i in range(1, N_total + 1):
        IJKi = essym.Num2IJK(i, N, Ord)
        IJK_list.append(IJKi)
    ijk_voigt = cij_voigt(IJK_list)
    (ijk_voigt_uniq, ijk_uniq_index) = np.unique(ijk_voigt, return_inverse=True)
    return (ijk_voigt_uniq, ijk_uniq_index)

def gen_coefmat_2d(symopt, Ord, N=2):
    #symopt is the symmetry operator, 2X2
    #Ord is the order (True order, not Voigt's notation)
    #
    #Solution: turn multi-loop to a single loop by 3-based number
    #    Turn Cijkl = aim*ajn*ako*alp*Cmnop
    #    To be C = Aij*C, here C is a column vector => Bij = 0 (Bij = Aij - I)
    #N = 2  # 2-based
    N_total = N ** Ord
    A = np.ones((N_total, N_total))
    for i in range(1, N_total + 1):
        # in Num2IJK (in symmetry), the first number is 1, and the maximum is N**Ord
        IJKi = essym.Num2IJK(i, N, Ord)
        #ijk_voigt = cij_voigt(IJKi)
        for j in range(1, N_total + 1):
            IJKj = essym.Num2IJK(j, N, Ord)
            for k in range(0, len(IJKj)):
                aij = symopt[IJKi[k]-1, IJKj[k]-1]
                A[i-1, j-1] = A[i-1, j-1] * aij
    B = A - np.eye(N_total)
    return B

def gen_coefmat_point_2d(symkey, Ord, N=2, tor=1.e-5):
    '''
    Generate the coefficient matrix for specific point group in 2d materials

    Parameter
    ---------
        symkey: str
            The point group. one of {'2', 'mm2', '4mm', '6mm'}
        Ord: ind
            The order of elastic constants in tensor order. For SOECs, it's 4
        N: ind
            The dimension of the materials, N=2 for 2d materials
        tor: float
            The tolerance for 0
    Return
    ------
        ijk_voigt_uniq: list
            The all unique cij. e.g. [111, 112, 116, 122, 126, 166, 222, 226, 266, 666]
        B1: list
            The elements that equal to ZERO
        B2: 2d list
            The expression for dependent elastic constants
    '''
    (ijk_voigt_uniq, ijk_uniq_index) = get_cij_uniq_from_ord(Ord=Ord, N=N)
    n = len(ijk_uniq_index)
    n_uniq = len(ijk_voigt_uniq)
    optnum = symoptgendb_2d(symkey)
    #B = []
    symopts = get_full_symopt_2d(symkey, tor=tor)
    i = 0
    #Get all the equations in tensor format
    for symopt in symopts:
        B0 = gen_coefmat_2d(np.array(symopt), Ord, N=N)
        if i == 0:
            B = B0
        else:
            B = np.vstack((B, B0))
        i = i + 1
    #Convert tensor to Voigt (sum for the same voigt)
    B_voigt = np.zeros((len(B), n_uniq))
    for i in range(n):
        i_voigt = ijk_uniq_index[i]
        B_voigt[:, i_voigt] = B_voigt[:, i_voigt] + B[:, i]

    #Get the number of independent equations by rank
    B_voigt[np.abs(B_voigt) < tor] = 0
    rnkB = np.linalg.matrix_rank(B_voigt)

    #Split to two parts: B1 (only one elements, =0), B2 (relations)
    (B1, B2) = split_B_by_N_NonZero(B_voigt, Ntor=1)
    B1 = get_fullrank(B1)
    Zero_ind = []
    B1_out = []
    for j, B1i in enumerate(B1):
        for i, b1i in enumerate(B1i):
            if b1i != 0:
                Zero_ind.append(i)
                B1_out.append(ijk_voigt_uniq[i])
                break
    rnkB1 = len(B1_out)
    #Set B2[:, B1] = 0
    for i, B2i in enumerate(B2):
        for j in Zero_ind:
            B2[i][j] = 0.
    '''    
    (B21, B22) = split_B_by_N_NonZero(B2, Ntor=3)
    B21 = get_fullrank(B21)
    rnkB21 = np.linalg.matrix_rank(B21)

    if rnkB21 == 0:
        B2 = get_fullrank(B22)
    elif rnkB21 == rnkB - rnkB1:
        B2 = B21
    else:
        #TODO
        B2 = get_fullrank(B2)
    '''
    #Normalize B2
    B2_out = []
    B2 = get_fullrank(B2)
    for B2i in B2:
        B2_out.append(norm_equation(B2i, tor=1e-5))

    return (ijk_voigt_uniq, B1_out, np.array(B2_out))

def simplify_B2(ijk_voigt_uniq, B1, B2):
    '''
    Simplify the expression
    
    '''
    (ijk_ind, ijk_eq, ijk_norm, ijk_shear) = determine_indep(ijk_voigt_uniq, B1, B2)
    ijk_noInd = np.setdiff1d(ijk_voigt_uniq, B1 + ijk_ind)
    n_noInd = len(ijk_noInd)

    coefmat_noInd = get_B_by_column(ijk_voigt_uniq, ijk_noInd, B2)
    B2_tmp = np.hstack((coefmat_noInd, B2))
    
    #print(B2_tmp)
    tor = 1e-5
    for j in range(n_noInd):
        for i in range(j, n_noInd):
            if i == j:
                k = j + 1
                while abs(B2_tmp[j][j]) < tor:
                    #If B2_tmp[j][j] == 0
                    B2_tmp[[j, k], :] = B2_tmp[[k, j], :]
                    k = k + 1
                B2_tmp[j] = B2_tmp[j] / B2_tmp[j][j]
            else:
                B2_tmp[i] = B2_tmp[i] - B2_tmp[j] * B2_tmp[i][j]
    
    for j in range(n_noInd):
        #n_noInd - j is the real row to be substract
        for i in range(j+1, n_noInd):
            i_real = n_noInd - i - 1
            j_real = n_noInd - j - 1
            B2_tmp[i_real] = B2_tmp[i_real] - B2_tmp[j_real] * B2_tmp[i_real][j_real]

    if len(B2_tmp) == 0:
        B2 = []
    else:
        B2 = B2_tmp[:, n_noInd:]
    
        B2[np.abs(B2)<1e-5] = 0
    
    return (ijk_ind, B2)

def get_B_by_column(ijk_voigt_uniq, ijk_column, B):
    '''
    Get the sub set of B by ijk_column
    '''
    coefmat = []
    for i, ijk_i in enumerate(ijk_voigt_uniq):
        if ijk_i in ijk_column:
            coefmat.append(B[:, i])
    coefmat = np.array(coefmat).T
    return coefmat

def determine_indep(ijk_voigt_uniq, B1, B2):
    '''
    Determining the independent elements.
    NOTE:
        The principle for selecting independent elastic constants
            1. The all-equal elements (e.g. 111, 222) should be the first.
            2. The no-shear components (e.g. 112, 123)
            3. The shear related (e.g. 116, 456)
    Parameter
    ---------
        ijk_voigt_uniq: list
            The all unique cij. e.g. [111, 112, 116, 122, 126, 166, 222, 226, 266, 666]
        B1: list
            The elements that equal to ZERO
        B2: 2d list
            The expression for dependent elastic constants
    Return
    ------
        ijk_ind: list
            The independent ijk
        ijk_eq: list
            The equal ijk
        ijk_norm: list
            The normal ijk (no shear-related component)
        ijk_shear: list
            The shere-related ijk
    '''
    def ijk2ten(ijk):
        ijk = str(ijk)
        ten = [int(ijk_i) for ijk_i in ijk]
        return ten
    def isequal(ten, tor=1e-5):
        delta = np.array(ten) - ten[0]
        if sum(np.abs(delta)) < tor:
            return True
        else:
            return False
    def isshear(ten):
        for teni in ten:
            if teni in [4, 5, 6]:
                return True
        else:
            return False
    
    rnkB2 = np.linalg.matrix_rank(B2)
    n_ind = len(ijk_voigt_uniq) - len(B1) - rnkB2
    ijk_eq, ijk_shear, ijk_norm, ijk_shear_eq = [], [], [], []
    for ijk in ijk_voigt_uniq:
        if ijk in B1:
            #equal to 0
            pass
        else:
            ten = ijk2ten(ijk)
            if isshear(ten):
                if isequal(ten):
                    ijk_shear_eq.append(ijk)
                else:
                    ijk_shear.append(ijk)
            elif isequal(ten):
                ijk_eq.append(ijk)
            else:
                ijk_norm.append(ijk)
    ijk_all = ijk_eq + ijk_norm + ijk_shear_eq + ijk_shear
    ijk_ind_list = list(combinations(ijk_all, n_ind))
    for ijk_columni in ijk_ind_list:
        ijk_column = np.setdiff1d(ijk_voigt_uniq, B1 + list(ijk_columni))
        coefmat = get_B_by_column(ijk_voigt_uniq, ijk_column, B2)
        if np.linalg.matrix_rank(coefmat) == rnkB2:
            ijk_ind = list(ijk_columni)
            break
    #ijk_ind = ijk_all[0:n_ind]
    #print(ijk_all)
    return (ijk_ind, ijk_eq, ijk_norm, ijk_shear)

def norm_equation(eq_list, tor=1e-5):
    '''
    Normalize the equation according to the minimum (abs of non-zero elements)
    Parameter
    ---------
        eq_list: list
            The coefficient list of equations
        tor: float
            The tolerance for zero elements
    Return
    ------
        eq_list: list
            The nomalized equation list
    '''
    eq_list_tmp = []
    #eq_list = 
    eq_abs_noZero = [abs(eqi) for eqi in eq_list if abs(eqi) > tor]
    eq_min = min(eq_abs_noZero)
    for eqi in eq_list:
        eq_list_tmp.append(0 if abs(eqi) < tor else eqi)
    eq_list_tmp = eq_list_tmp / eq_min    
    return eq_list_tmp

def count_non_zeros(Bi, tor=1e-5):
    b = [1 for bi in Bi if abs(bi) > tor]
    return sum(b)

def split_B_by_N_NonZero(B, Ntor=6, tor=1e-5):
    '''
    Split B according to the number of non-zeros elements (Ntor)
    '''
    B1, B2 = [], []
    for Bi in B:
        n = count_non_zeros(Bi)
        if n == 0:
            pass
        elif n > Ntor:
            B2.append(Bi)
        else:
            B1.append(Bi)
    return (B1, B2)

def get_fullrank(B):
    # Remove some vectors which is linear dependent
    #(B1, B2) = split_B_by_N_NonZero(B)
    #rnkB = np.linalg.matrix_rank(B)
    B = np.array(B)
    B[np.abs(B)<1e-5] = 0
    rnk0 = np.linalg.matrix_rank(B)
    if rnk0 == 1:
        return np.reshape(B[0], (1, -1))
    for i in range(0, len(B)):
        rnk0 = np.linalg.matrix_rank(B)
        #print(rnk0)
        B0 = B[0]
        B = np.delete(B, 0, 0)
        rnk1 = np.linalg.matrix_rank(B)
        if rnk0 != rnk1:
            B = np.append(B, np.reshape(B0, (1, -1)), axis=0)
    return B

def get_full_symopt_2d(symkey, tor=1e-5):
    '''
    Generate the full symmetry operation by specific the symmetry key
        It generated by the following algrithm:
            1. Get the generate elements from the db
            2. For each elements, get the full symmetries (symopt x symopt x symopt ... = e)
            3. Remove the redundent symopt
            4. Combine each other (Order sensitive)
    Parameter
    ---------
        symkey: str
            The symmetry key for 2d lattice, one of {'2', 'mm2', '4mm', '6mm'}
        tor: float
            The tolerance of redundent
    Return
    ------
        symopt_full: list, nx2x2
            The full symmetry operations for specific symkey. nx2x2
    '''
    optnum = symoptgendb_2d(symkey)
    symopts, symopt_full = [], []
    for i in range(0, len(optnum)):
        optmat = symoptdb_2d(optnum[i])
        symopts = symopts + get_full_symopt_single_2d(optmat)
    symopts = rem_redun_symopt_2d(symopts, tor=tor)
    for symopti in symopts:
        for symoptj in symopts:
            symopt_full.append(np.dot(symopti, symoptj).tolist())
    symopt_full = rem_redun_symopt_2d(symopt_full, tor=tor)
    return symopt_full

def rem_redun_symopt_2d(symopts, tor=1e-5):
    '''
    Remove the redundent elements in the symmetry operations
        Compare the sum of the absolute differents of two symmetry operations
    Parameter
    ---------
        symopts: list, nx2x2
            The symmetry operations with/without redundent
        tor: float
            The tolerance of redundent
    Return
    ------
        symopt_r: list, nx2x2
            The symmetry operations without redundent. 
    '''
    n = len(symopts)
    symopt_r = []
    for i in range(0, n):
        symopti = np.array(symopts[i])
        for j in range(i+1, n):
            symoptj = np.array(symopts[j])
            deltaij = cal_delta_2d(symopti, symoptj)
            if deltaij < tor:
                break
        else:
            symopt_r.append(symopti.tolist())
    return symopt_r        

def cal_delta_2d(symopt, symopt2=np.eye(2)):
    '''
    Calculate the difference of two symmetry operations
        It calculate the sum of the absolute differents of two symmetry operations
        If symmat2 is not specific, then compare with the np.eye(2)
    Parameter
    ---------
        symopt: 2x2 list
            The symmetry operations for 2d materials
        symopt2: 2x2 list
            The symmetry operations for 2d materials
    Return
    ------
        delta: float
            The difference of two symmtry operations.
    '''
    symopt_diff = symopt - symopt2
    return np.sum(np.abs(symopt_diff))

def get_full_symopt_single_2d(symopt, tor=1e-5):
    '''
    Get all the symmetry operations for a specif symmetry operation generator (x, x^2, x^3, x^4, ..., e)
    Parameter
    ---------
        symopt: list, nx2x2
            The symmetry operations with/without redundent
        tor: float
            The tolerance of redundent
    Return
    ------
        symopts: list, nx2x2
            The symmetry operations for a specific symmetry operation generator. 
    '''
    symopt = np.array(symopt)
    symopt0 = copy.deepcopy(symopt)
    symopts = []
    while cal_delta_2d(symopt) > tor:
        symopts.append(symopt.tolist())
        symopt = np.dot(symopt, symopt0)
    symopts.append(np.eye(2).tolist())
    return symopts

def get_relationship(symkey, Ord):
    (ijk_voigt_uniq, B1, B2) = gen_coefmat_point_2d(symkey, Ord)
    (ijk_ind, B2) = simplify_B2(ijk_voigt_uniq, B1, B2)

    print('The point group is: {}\n'.format(symkey))
    print('The All IJK List are:\n{}\n'.format(ijk_voigt_uniq))
    print('The Zero IJK are: \n{}\n'.format(B1))
    print('The Independent IJK are: \n{}\n'.format(ijk_ind))
    print('The Relationship between Dependent IJK are: \n{}\n{}\n'.format(ijk_voigt_uniq, B2))
    return (ijk_voigt_uniq, ijk_ind, B1, B2)

'''
symkey = '6mm'
symopts = get_full_symopt_2d(symkey)
print(symopts)
print(len(symopts))

Ord = 4
symopt = np.array([[1., 0.], [0., -1.]])
B = gen_coefmat_2d(symopt, Ord)
print(B)

Ord = 8
symkey = '6mm'
(ijk_voigt_uniq, B1, B2) = gen_coefmat_point_2d(symkey, Ord)
print(ijk_voigt_uniq)
print(B1)

(ijk_ind, B2) = simplify_B2(ijk_voigt_uniq, B1, B2)
print(B2)
print(ijk_ind)
'''

Ord = 12
symkey = '6mm'
get_relationship(symkey, Ord)