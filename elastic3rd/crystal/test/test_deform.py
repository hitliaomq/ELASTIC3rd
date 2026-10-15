#!python

import elastic3rd.crystal.deform as crydef
import numpy as np

StrainVerctor = np.array([1.0, 0.0, 1.0, 2.0, 2.0, 0.0])
StrainMatrix = np.array([[1.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]])

def test_vec2matrix_matrix2vec():
    StrainMatrix1 = crydef.vec2matrix(StrainVerctor)
    assert(StrainMatrix == StrainMatrix1).all()

    StrainVerctor1 = crydef.matrix2vec(StrainMatrix)
    assert(StrainVerctor == StrainVerctor1).all()

def test_is_strain_matrix():
    IsStrainMatrix = crydef.is_strain_matrix(StrainMatrix)
    assert(IsStrainMatrix)

    StrainMatrix_Wrong1 = np.array([[1, 2], [2, 3]])
    NotStrainMatrix1 = crydef.is_strain_matrix(StrainMatrix_Wrong1)
    assert(not NotStrainMatrix1)
    StrainMatrix_Wrong2 = np.array([[1, 2, 3], [2, 2, 3], [3, 2, 2]])
    NotStrainMatrix2 = crydef.is_strain_matrix(StrainMatrix_Wrong2)
    assert(not NotStrainMatrix2)

def test_strain2deformgrad():
    ActualStrainMatrix = 0.03 * StrainMatrix
    deformgrad1 = crydef.strain2deformgrad(ActualStrainMatrix)
    deformgrad2 = crydef.strain2deformgrad2(ActualStrainMatrix)
    assert(np.allclose(deformgrad1, deformgrad2, atol=1e-3))

def test_eigv():
    (Dv, V) = crydef.eigv(StrainMatrix)
    assert(crydef.is_strain_matrix(Dv))
    assert(crydef.is_strain_matrix(V))
