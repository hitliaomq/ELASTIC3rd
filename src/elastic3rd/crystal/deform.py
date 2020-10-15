#/usr/bin/python
#

import numpy as np

def deform_mode(strain, flag):
    '''
    Generate strain mode for cubic and trigonal system.
    (Note: This is the initial version of Elastic3rd, not used in current version)
             ref:PHYSICAL REVIEW B 75, 094105 (2007)
    Parameters
    ----------
        strain: float
            The amplitude of the strain
        flag: int
            The index of the strain mode.
    Return
    ------
        STRAIN: np.ndarray
            The list of strian, Voigt notation
    '''
    strain = float(strain)
    if flag == 1:
        STRAIN = np.array([strain, 0, 0, 0, 0, 0])
    elif flag == 2:
        STRAIN = np.array([strain, strain, 0, 0, 0, 0])
    elif flag == 3:
        STRAIN = np.array([strain, strain, strain, 0, 0, 0])
    elif flag == 4:
        STRAIN = np.array([strain, 0, 0, 2. * strain, 0, 0])
    elif flag == 5:
        STRAIN = np.array([strain, 0, 0, 0, 0, 2. * strain])
    elif flag == 6:
        STRAIN = np.array([0, 0, 0, 2. * strain, 2. * strain, 2. * strain])
    elif flag == 7:
        STRAIN = np.array([0, strain, 0, 0, 0, 0])
    elif flag == 8:
        STRAIN = np.array([0, 0, strain, 0, 0, 0])
    elif flag == 9:
        STRAIN = np.array([0, strain, strain, 0, 0, 0])
    elif flag == 10:
        STRAIN = np.array([0, 0, 0, 2. * strain, 0, 0])
    elif flag == 11:
        STRAIN = np.array([0, 0, strain, 0, 0, 2. * strain])
    elif flag == 12:
        STRAIN = np.array([0, strain, 0, 2. * strain, 0, 0])
    elif flag == 13:
        STRAIN = np.array([0, strain, 0, 0, 2. * strain, 0])
    elif flag == 14:
        STRAIN = np.array([0, 0, strain, 0, 2. * strain, 0])
    elif flag == 15:
        STRAIN = np.array([strain, 0, strain, 2. * strain, 0, 0])
    elif flag == 16:
        STRAIN = np.array([strain, strain, 0, 0, 2. * strain, 0])

    return STRAIN

def vec2matrix(StrainVerctor):
    '''
    Transform the strain list from verctor to matrix (All in Voigt's notation)
    Parameters
    ----------
        StrainVerctor: list/tuple
            The strain in verctor format, e.g. [1.0, 0.0, 1.0, 2.0, 2.0, 0.0]
    Return
    ------
        StrainMatrix: 2D np.ndarray
            The strain in matrix format, e.g.
            [[1.0, 0.0, 1.0],
             [0.0, 0.0, 1.0],
             [1.0, 1.0, 1.0]]
    '''
    (e1, e2, e3, e4, e5, e6) = StrainVerctor
    e4 = e4/2.0
    e5 = e5/2.0
    e6 = e6/2.0
    StrainMatrix = np.array([[e1, e6, e5], [e6, e2, e4], [e5, e4, e3]])
    return StrainMatrix

def matrix2vec(StrainMatrix):
    '''
    Transform the strain list from matrix to vector (All in Voigt's notation)
    Parameters
    ------
        StrainMatrix: 2D np.ndarray
            The strain in matrix format, e.g.
            [[1.0, 0.0, 1.0],
             [0.0, 0.0, 1.0],
             [1.0, 1.0, 1.0]]
    Return
    ----------
        StrainVerctor: list/tuple
            The strain in verctor format, e.g. [1.0, 0.0, 1.0, 2.0, 2.0, 0.0]
    '''
    if not is_strain_matrix(StrainMatrix):
        StrainVerctor = np.array([0, 0, 0, 0, 0, 0])
        print("Warning: The format of input strain matrix is not correct. And it is set as zeros")
        return StrainVerctor
    e1 = StrainMatrix[0][0]
    e2 = StrainMatrix[1][1]
    e3 = StrainMatrix[2][2]
    e4 = 2. * StrainMatrix[1][2]
    e5 = 2. * StrainMatrix[0][2]
    e6 = 2. * StrainMatrix[0][1]
    StrainVerctor = np.array([e1, e2, e3, e4, e5, e6])
    return StrainVerctor

def is_strain_matrix(StrainMatrix):
    '''
    Judge if the given matrix is a strain matrix or not, accordint to two condition
        1. size = 3x3
        2. symmetrical
    Parameters
    ----------
        StrainMatrix: 2D np.ndarray
            The strain in matrix format, e.g.
            [[1.0, 0.0, 1.0],
             [0.0, 0.0, 1.0],
             [1.0, 1.0, 1.0]]
        tolerance: float
            The tolerance for symmetrical
    Return
    ------
        flag: bool
            If StrainMatrix is correct, return True, ortherwise return False
    '''
    (m, n) = StrainMatrix.shape
    if not (m == 3 and n == 3):
        print("Error: The size of input strain matrix is not 3 by 3!")
        return False
    else:
        for i in range(1, 3):
            for j in range(0, i):
                if StrainMatrix[i][j] != StrainMatrix[j][i]:
                    print("Warning: The input strain matrix is not symmtry!")
                    return False
    return True

def strain2deformgrad(StrainMatrix):
    '''
    Calculate the deformation gradient according to Lagrangian strain 
        (Note: only work for symmetrical strain matrix)
    Parameters
    ----------
        StrainMatrix: np.ndarray
            The Lagrangian strain in matrix format
    Return
    ------
        DeformGrad: float
            The deformation gradient
    '''
    Y = 2 * StrainMatrix + 1 * np.eye(3)
    (D, V) = eigv(Y)
    VT = V.T
    DeformGrad = V.dot(np.sqrt(D)).dot(VT)
    return DeformGrad

def strain2deformgrad2(StrMat):
    '''
    Calculate the deformation gradient according to Lagrangian strain
        Ref: Eq.6 in Acta Materialia 142 (2018)  [10.1016/j.actamat.2017.09.041]
    Parameters
    ----------
        StrainMatrix: np.ndarray
            The Lagrangian strain
    Return
    ------
        DeformGrad: float
            The deformation gradient
    '''
    DeformGrad = np.zeros((3, 3))
    for i in range(0, 3):
        for j in range(0, 3):
            if i == j:
                deltaij = 1
            else:
                deltaij = 0
            DeformGrad[i, j] = deltaij + StrMat[i, j]
            for k in range(0, 3):
                DeformGrad[i, j] = DeformGrad[i, j] - 0.5*StrMat[k, i]*StrMat[k, j]
                for l in range(0, 3):
                    DeformGrad[i, j] = DeformGrad[i, j] + 0.5*StrMat[k, i]*StrMat[l, k]*StrMat[l, j]
    return DeformGrad

def eigv(a):
    '''
    Calculate the eigenvalues and eigenvectors for symmetrical matrix
    Parameter
    ---------
        a: 2D np.ndarray
            The matrix to calculate the eig value and eig vector
    Returns
    -------
        Dv: 2D np.ndarray
            The eigenvalues matrix. The eig value loacted along the diag
        V: 2D np.ndarray
            The eigenvectors matrix
    '''
    (D, V) = np.linalg.eigh(a)
    #print D
    #print V
    Dv = np.zeros((3, 3))
    for i in range(0, 3):
        Dv[i, i] = D[i]
    return (Dv, V)
