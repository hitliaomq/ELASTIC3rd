#!python
#
#symopt: symmetry operater
#
#NOTE: This script is not used now. It's for n-th order elastic constants usage
#
import elastic3rd.symmetry.symmetry as essym
import numpy as np
import copy
from scipy.linalg import lu
import math

def symopt(Dim = 3, sym = '1', axis = 'x'):
    # Dim is the dimension
    # sym is the independent symmetry operation, here it equals to 1, i, m, 2, 3, 4, 6 and -4
    # axis is a verctor
    #   sym = 1 or i, no meaning
    #   sym = m, axis is the nomal verctor of plane
    #   sym = 2, 3, 4, 6 or -4, axis is the axis of rotation
    if Dim == 2:
        if sym == '1':
            TransfMatrix = np.array([[1, 0], [0, 1]])
        elif sym == 'i' or sym == 'I':
            TransfMatrix = np.array([[-1, 0], [0, -1]])
        elif sym == '2' or sym == '3' or sym == '4' or sym == '6':
            #ang = 0
            if sym == '2':
                ang = 180.
            elif sym == '3':
                ang = 120.
            elif sym == '4':
                ang = 90.
            elif sym == '6':
                ang = 60.
            ang = np.radians(ang)
            TransfMatrix = np.array([[math.cos(ang), -math.sin(ang)], [math.sin(ang), math.cos(ang)]])
        elif sym == 'm':
            if axis == 'x':
                TransfMatrix = np.array([[-1, 0], [0, 1]])
            elif axis == 'y':
                TransfMatrix = np.array([[1, 0], [0, -1]])
            else:
                raise Exception("No such mirror plane")
        else:
            raise Exception("No such symmetry in sym")
    elif Dim == 3:
        if sym == '1':
            TransfMatrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        elif sym == 'i' or sym == 'I':
            TransfMatrix = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        elif sym == '2' or sym == '3' or sym == '4' or sym == '6' or sym == '-4':
            #ang = 0
            if sym == '2':
                ang = 180.
            elif sym == '3':
                ang = 120.
            elif sym == '4' or sym == '-4':
                ang = 90.
            elif sym == '6':
                ang = 60.
            ang = np.radians(ang)
            axis = get_axis(axis)
            TransfMatrix = RotationMatrix(axis, ang)
            if sym == '-4':
                # add a inverse on 4 for -4
                TransfMatrix = -TransfMatrix
        elif sym == 'm':
            axis = get_axis(axis)
            if (axis == np.array([1, 0, 0])).any():
                TransfMatrix = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
            elif (axis == np.array([0, 1, 0])).any():
                TransfMatrix = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
            elif (axis == np.array([0, 0, 1])).any():
                TransfMatrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
            else:
                raise Exception("No such mirror plane")
        else:
            raise Exception("No such symmetry in sym")
    return TransfMatrix

def get_axis(axis):
    # convert axis to direction
    sqrt3 = np.sqrt(3.)/3.  # for normalizing the 111 direction
    if type(axis) is str:
        # for some special axis, eg : 'x', 'y', 'z'
        if axis.lower() == 'x' or axis == '1':
            axis = np.array([1, 0, 0])
        elif axis.lower() == 'y' or axis == '2':
            axis = np.array([0, 1, 0])
        elif axis.lower() == 'z' or axis == '3':
            axis = np.array([0, 0, 1])
        elif axis.lower() == 'p1':
            axis = sqrt3 * np.array([1, 1, 1])
        elif axis.lower() == 'p2':
            axis = sqrt3 * np.array([1, 1, -1])
        elif axis.lower() == 'p3':
            axis = sqrt3 * np.array([1, -1, 1])
        elif axis.lower() == 'p4':
            axis = sqrt3 * np.array([1, -1, -1])
    else:
        # it must be a verctor with 3 numerical elements
        #   then normalize the vector
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
    return axis

def RotationMatrix(m, ang):
    #axis: x-x axis; y-y axis; z-z axis; p-[111] axis
    #m: direction in angle cosine, eg: x axis [1, 0, 0]
    #ang: angle in degree
    #ang = math.radians(ang)    
    cosa = math.cos(ang)
    sina = math.sin(ang)
    cos1 = 1 - cosa
    m1 = m[0]
    m2 = m[1]
    m3 = m[2]
    R_Maxtix = np.array([[cosa+m1*m1*cos1, m1*m2*cos1-m3*sina, m1*m3*cos1+m2*sina], 
        [m2*m1*cos1+m3*sina, cosa+m2*m2*cos1, m2*m3*cos1-m1*sina], 
        [m3*m1*cos1-m2*sina, m3*m2*cos1+m1*sina, cosa+m3*m3*cos1]])
    return R_Maxtix

def symopt_point(pointsym):
    #ROBERT E. NEWNHAM. Properties of Materials :Anisotropy, Symmetry, Structure, P24
    #   Oxford University Press, 2005
    sym = pointsym.lower()
    if sym == "1":
        pass

def symoptgendb(symkey):
    symoptdb_var = {'1' : [1], '-1' : [2], '2' : [3], 'm': [4], '2/m' : [3, 4], '222' : [5, 3], 'mm2' : [6, 3], 
                'mmm' : [6, 7, 4], '4' : [8], '-4' : [9], '422' : [8, 5], '4/m' : [8, 4], '4mm' : [8, 6], 
                '-42m' : [9, 5], '4/mmm' : [8, 4, 6], '3' : [10], '-3' : [11], '32' : [10, 5], '3m' : [10, 6], 
                '-3m' : [11, 6], '6' : [12], '-6' : [13], '-6m2' : [13, 6], '622' : [12, 5], '6/m' : [12, 4], 
                '6mm' : [12, 6], '6/mmm' : [12, 4, 6], '23' : [3, 14], '432' : [8, 14], 'm-3' : [4, 15], 
                '-43m' : [9, 14],'m-3m' : [8, 15, 16]}
    optnum = symoptdb_var[symkey]
    return optnum

def symoptdb(NO):
    #Ref: Vainsthein B.K. (1994) Fundamentals of the Theory of Symmetry. In: Fundamentals of Crystals. Springer, Berlin, Heidelberg P98-99
    #      Correct the sqrt(3)/2.0 (Original is sqrt(3/2))
    #Note : this is only the basic transformation matrix, it need to be expand
    sqrt3 = math.sqrt(3.)/2.0
    symdb = {
        1 : np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 
        2 : np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]),
        3 : np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]), 
        4 : np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),
        5 : np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), 
        6 : np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        7 : np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]]), 
        8 : np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
        9 : np.array([[0, -1, 0], [1, 0, 0], [0, 0, -1]]), 
        10 : np.array([[-0.5, sqrt3, 0], [-sqrt3, -0.5, 0], [0, 0, 1]]),
        11 : np.array([[0.5, -sqrt3, 0], [sqrt3, 0.5, 0], [0, 0, -1]]), 
        12 : np.array([[0.5, sqrt3, 0], [-sqrt3, 0.5, 0], [0, 0, 1]]), 
        13 : np.array([[-0.5, -sqrt3, 0], [sqrt3, -0.5, 0], [0, 0, -1]]), 
        14 : np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
        15 : np.array([[0, -1, 0], [0, 0, -1], [-1, 0, 0]]), 
        16 : np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    }
    optmat = symdb[NO]
    return optmat

def symoptdb2(NO):
    #NOT USED
    #ROBERT E. NEWNHAM. Properties of Materials :Anisotropy, Symmetry, Structure, P24
    #   Oxford University Press, 2005
    #   the NO.16 is missing
    sqrt3 = math.sqrt(3.)/2.
    symdb = {
        '1' : np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        '2' : np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]),
        '3' : np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),
        '4' : np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]),
        '5' : np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        '6' : np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]]),
        '7' : np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),
        '8' : np.array([[-0.5, sqrt3, 0], [-sqrt3, -0.5, 0], [0, 0, 1]]),
        '9' : np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
        '10': np.array([[0.5, -sqrt3, 0], [sqrt3, 0.5, 0], [[0, 0, -1]]]),
        '11': np.array([[0, -1, 0], [0, 0, -1], [-1, 0, 0]]),
        '12': np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
        '13': np.array([[0, -1, 0], [1, 0, 0], [0, 0, -1]]),
        '14': np.array([[0.5, sqrt3, 0], [-sqrt3, 0.5, 0], [0, 0, 1]]),
        '15': np.array([[-0.5, -sqrt3, 0], [sqrt3, -0.5, 0], [0, 0, -1]]),
        '16': np.array([[1,1,1,1]])
    }

def point2laue(point):
    laue = {'-1' : ['1', '-1'], '2/m' : ['2', 'm', '2/m'], 'mmm' : ['222', 'mm2', 'mmm'], '4/m' : ['4', '-4', '4/m'], 
    '4/mmm' : ['422', '4mm', '-42m', '4/mmm'], '-3' : ['3', '-3'], '-3m' : ['32', '3m', '-3m'], '6/m' : ['6', '-6', '6/m'],
    '6/mmm' : ['622', '6mm', '-6m2', '6/mmm'], 'm-3' : ['23', 'm-3'], 'm-3m' : ['432', '-43m', 'm-3m']}
    for key in laue:
        valuei = laue[key]
        if point in valuei:
            print(key)
            return key
    print('This is not correct point group.')

def laue2crystaltype(laue):
    lauecrytyp = {'-1' : 'n', '2/m' : 'm', 'mmm' : 'o', '4/m' : 't2', '4/mmm' : 't1', '-3' : 'r2',
    '-3m' : 'r1', '6/m' : 'h2', '6/mmm' : 'h1', 'm-3' : 'c2', 'm-3m' : 'c1'}
    CrystalType = lauecrytyp[laue]
    return CrystalType

def crystaltype2laue(CrystalType):
    crytyplaue = {'n' : '-1', 'm' : '2/m', 'o' : 'mmm', 't2' : '4/m', 't1' : '4/mmm', 'r2' : '-3',
    'r1' : '-3m', 'h2' : '6/m', 'h1' : '6/mmm', 'c2' : 'm-3', 'c1' : 'm-3m'}
    laue = crytyplaue[CrystalType]
    return laue

def get_equation(symmat, Ord):
    orishape = symmat.shape
    for i in range(0, Ord/2):
        symmat = np.outer(symmat, symmat)
    symmat = symmat.reshape(Ord * orishape)
    print(symmat.shape)
    axes = (range(0, 2*Ord, 2), range(0, Ord))
    print(axes)
    a = np.tensordot(symmat, np.ones(Ord*orishape), axes)
    print(a.shape)
    print(a)
    return a

def gen_coefmat(optmat, Ord):
    #symopt is the symmetry operator, 3x3
    #Ord is the order (True order, not Voigt's notation)
    #
    #Solution: turn multi-loop to a single loop by 3-based number
    #    Turn Cijkl = aim*ajn*ako*alp*Cmnop
    #    To be C = Aij*C, here C is a column vector => Bij = 0 (Bij = Aij - I)
    N = 3  # 3-based
    N_total = N ** Ord
    A = np.ones((N_total, N_total))
    for i in range(1, N_total + 1):
        # in Num2IJK (in symmetry), the first number is 1, and the maximum is N**Ord
        IJKi = essym.Num2IJK(i, N, Ord)
        for j in range(1, N_total + 1):
            IJKj = essym.Num2IJK(j, N, Ord)
            for k in range(0, len(IJKj)):
                aij = optmat[IJKi[k]-1, IJKj[k]-1]
                A[i-1, j-1] = A[i-1, j-1] * aij
    #print(A)
    B = A - np.eye(N_total)
    return B

def gen_coefmat_point(symkey, Ord):
    optnum = symoptgendb(symkey)
    #B = []
    for i in range(0, len(optnum)):
        optmat = symoptdb(optnum[i])
        B0 = gen_coefmat(optmat, Ord)
        if i == 0:
            B = B0
        else:
            B = np.vstack((B, B0))
    #B = rem_redun(B)
    #B=list(set([tuple(t) for t in B]))
    B = get_fullrank(B)
    return B

def get_full_symopt(symkey):
    optnum = symoptgendb(symkey)
    symopts = []
    for i in range(0, len(optnum)):
        optmat = symoptdb(optnum[i])
        symopts = symopts + get_full_symopt_single(optmat)
    return symopts

def rem_redun_symopt(symopts, tor=1e-5):
    n = len(symopts)
    symopt_r = []
    for i in range(0, n):
        symopti = np.array(symopts[i])
        flag = 1
        for j in range(i+1, n):
            symoptj = np.array(symopts[j])
            deltaij = cal_delta(symopti, symoptj)
            if deltaij < tor:
                flag = 0
                break
        if flag:
            symopt_r.append(symopti.tolist())
    return symopt_r        

def cal_delta(symmat, symmat2=np.eye(3)):
    symmat = symmat - symmat2
    return np.sum(np.abs(symmat))

def get_full_symopt_single(symmat, tor=1e-5):
    symmat = np.array(symmat)
    symmat0 = copy.deepcopy(symmat)
    symmats = []
    while cal_delta(symmat) > tor:
        symmats.append(symmat.tolist())
        symmat = np.dot(symmat, symmat0)
    symmats.append(np.eye(3).tolist())
    return symmats

def get_fullrank(B):
    # Remove some vectors which is linear dependent
    B = np.array(B)
    for i in range(0, len(B)):
        rnk0 = np.linalg.matrix_rank(B)
        #print(rnk0)
        B0 = B[0]
        B = np.delete(B, 0, 0)
        rnk1 = np.linalg.matrix_rank(B)
        if rnk0 != rnk1:
            B = np.append(B, np.reshape(B0, (1, -1)), axis=0)
    return B


def rem_redun(B):
    # remove the redundent equations
    # B=list(set([tuple(t) for t in B]))  # This is more efficient
    rem_lst = []
    #red_lst = []
    for i in range(1, len(B)):
        eqi = B[i]
        flag = 0  # to judge if it's redundent
        for j in range(0, i):
            eqj = B[j]
            for k in range(0, len(B[0])):
                #compare eqi with eqj for each elements
                if eqi[k] != eqj[k]:
                    flag = 0
                    break
                if k == len(B[0])-1:
                    flag = 1
            if flag == 1:
                break
        if flag == 1:
            rem_lst.append(i)
            #red_lst.append(j)
    return B

def simplfy_B(B):
    #Ideal: Take each equation as a_last = ....
    #         If some new relations are determined, update previous corresponding data
    #         What to do if current one is the same with previous one or conflict with previous ?
    #P, L, U =  lu(B)
    CoefInd = [i+1 for i in range(0, len(B[0]))]
    IndInd = [i+1 for i in range(0, len(B[0]))]
    CoefCoef = [1 for i in range(0, len(B[0]))]
    for i in range(0, len(B)):
        #in each equation
        eqi = B[i]
        Num_eq = []
        Ind_eq = []
        Coef_eq = []
        for j in range(0, len(eqi)):
            # Get the nonzeros elements, Num_eq is not the index, but the index + 1
            if eqi[j] != 0:
                Num_eq.append(CoefInd[j])
                Ind_eq.append(j)   #This is the ture index of the position of nonzeros elements
                Coef_eq.append(eqi[j])
        if Num_eq:
            # Not empty
            # 
            # find the index of the maximum
            max_ind = [i for i, j in enumerate(Num_eq) if j == max(Num_eq)]
            if len(max_ind) > 1:
                #TODO
                print("YES")
                max_ind = max_ind[0]
                pass
            else:
                max_ind = max_ind[0]
            #print(Ind_eq)
            Num_eq_pop = Ind_eq.pop(max_ind)  #take out the index of max value
            Num_eq.pop(max_ind) #delete the corresponding Num_eq
            Coef_eq_pop = Coef_eq.pop(max_ind)   #take out the max value
            for k in range(0, len(Num_eq)):
                #turn the equation to max_value = Coef_eq * Other_value
                Coef_eq[k] = float(Coef_eq[k])/(-float(Coef_eq_pop))
            if not Num_eq:
                CoefInd[Num_eq_pop] = 0
                CoefCoef[Num_eq_pop] = 0
                continue
            if len(Num_eq) > 1:
                CoefIndi = Num_eq
                CoefCoefi = Coef_eq
            else:
                #print(Num_eq)
                CoefIndi = Num_eq[0]
                CoefCoefi = Coef_eq[0]
            #CoefInd[Num_eq_pop] = CoefIndi
            #CoefCoef[Num_eq_pop] = CoefCoefi
            CoefInd_rep = CoefInd[Num_eq_pop]
            #CoefCoef_rep = CoefCoef[Num_eq_pop]
            for k in range(0, len(CoefInd)):
                #print(CoefInd)
                if CoefInd[k] == CoefInd_rep:
                    CoefInd[k] = CoefIndi
                    CoefCoef[k] = CoefCoef[k]*CoefCoefi
            #print(Num_eq)
            #print(Coef_eq)
    print(CoefInd)
    print(CoefCoef)


'''
symmat = np.array([[-0.5, math.sqrt(3.)/2., 0], [-math.sqrt(3.)/2., -0.5, 0], [0, 0, 1.]])
Ord = 2
get_equation(symmat, Ord)
a = range(0, Ord, 2)
print a
'''

'''
Dim = 3
sym = 'r'
flag = 'p'
ang = 120

TM = symopt(Dim, sym, flag, ang)
print(TM)
'''

Ord = 4

symkey = 'mm2'
B = gen_coefmat_point(symkey, Ord)

rnk = np.linalg.matrix_rank(B)
print(rnk)

#simplfy_B(B)

#np.savetxt('B.txt', B)
#simplfy_B(B)

symmat = symoptdb(16)
print(symmat)
symmats = get_full_symopt_single(symmat)
print(symmats)
for symmati in symmats:
    print(symmati)
print(len(symmats))


symkey = '1'
symopts = get_full_symopt(symkey)
print(symopts)
print(len(symopts))

symopt_r = rem_redun_symopt(symopts, tor=1e-5)
print(symopt_r)
print(len(symopt_r))
