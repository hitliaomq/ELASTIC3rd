#!python
#
#
from elastic3rd.symmetry.symdata import coef_crystal, coef_ind, get_unique
from elastic3rd.symmetry.symmetry import num2cijk


def indpen2full_list(IndepenEC, Ord=2, CrystalType="c1"):
    """
    Transfer the indpendent elastic constatnt to the full elastic constant list

    Parameters
    ----------
        IndepenEC: list
            Independent Elastic Constants
        Ord: int
            The order of the elastic constant
        CrystalType: str
            The crystal type
    Return
    ------
        FullEC: list
            The full list of the elastic constant 
    """
    FullEC = []
    ECind = coef_ind(CrystalType, Ord)
    PreCoef = coef_crystal(CrystalType, Ord)
    CoefUniq = get_unique(ECind)
    if len(CoefUniq) != len(IndepenEC):
        raise Exception("The number of indpendent elastic constant does not match the crystal type or order."
            + "The required input of indepenEC is " + str(len(CoefUniq)))
    for i in range(len(ECind)):
        if ECind[i] == 0:
            FullEC.append(0.)
        elif type(ECind[i]) is list:
            s = 0
            for j in range(len(ECind[i])):
                ind = CoefUniq.index(ECind[i][j])
                s += float(PreCoef[i][j] * float(IndepenEC[ind]))
            FullEC.append(s)
        else:
            ind = CoefUniq.index(ECind[i])
            FullEC.append(float(PreCoef[i] * float(IndepenEC[ind])))
    return FullEC

def indpen2full_dict(IndepenEC, Ord=2, CrystalType="c1"):
    """
    Transfer the indpendent elastic constatnt to the full elastic constant dict

    Parameters
    ----------
        IndepenEC: list
            Independent Elastic Constants
        Ord: int
            The order of the elastic constant
        CrystalType: str
            The crystal type
    Return
    ------
        FullEC_dict: dict
            The dict of the elastic constant, e.g. {"C11": 1054, "C12": 456, ...}
    """
    fullEC = indpen2full_list(IndepenEC, Ord, CrystalType)
    FullEC_dict = {}
    for i in range(len(fullEC)):
        FullEC_dict["C" + str(num2cijk(i+1, Ord))] = fullEC[i]
    return FullEC_dict 

def EC2matrix(fullEC, Ord=2, CrystalType="c1"):
    """
    Transfer the full_elastic_constant list to matrix

    Parameter
    ---------
        fullEC: list
            The list of full_elastic_constant
        Ord: int
            The order of the elastic constant
        CrystalType: str
            The crystal type
    Return
    ------
        EC_Matrix: np.ndarray
            The matrix of the elastic constant, with size 3x3x3x3(Ord=2) 3x3x3x3x3x3(Ord=3)
    """
    pass

def rotation():
    pass

'''
IndepenEC = [1111, 126, 562, 300, 100, 200, 124, 235, 453, 982]
Ord = 3
CrystalType = "h1"
FullEC = indpen2full_list(IndepenEC, Ord, CrystalType)
print(FullEC)

FullEC_dict = indpen2full_dict(IndepenEC, Ord, CrystalType)
print(FullEC_dict)
'''
