#!python
#

import numpy as np
import math
#import esfit as esfit
import elastic3rd.post.esfit as esfit
import elastic3rd.symmetry.symmetry as essym
import elastic3rd.symmetry.symmetry_stress as symstr
import elastic3rd.esutils as esutils
from elastic3rd.esutils import logger
import math

def get_cij(coef_fit, coef2, coef3, flag_se):
    '''
    Get the elastic constant by the fitting parameters
    Parameters
    ----------
        coef_fit: np.ndarray
            The fitting parameters
        coef2: np.ndarray
            The analytical coefficient matrix of SOECs
        coef3: np.ndarray
            The analytical coefficient matrix of TOECs
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
    Return
    ------
        C2: np.ndarray
            SOECs, solved by Least squares method
        C3: np.ndarray
            TOECs, if the number of coefficient equals to the number of independent equations, then solved by solving the equations
                   if the number of coefficient less than the number of independent equations, then solved by Least squares method
    '''
    flag_se = flag_se.lower()
    if flag_se == "e":
        C3 = np.linalg.solve(coef3, 6.0*coef_fit[:, 1])
        #C2 = np.linalg.solve(coef2[[0, 1, 6], [0, 1, 6]], 2.0*coef_fit[[0, 1, 6], 0])
        C2 = np.linalg.lstsq(coef2, 2.0*coef_fit[:, 0], rcond = None)[0]
        #print C2
        #print(residuals)
        np.savetxt("FittedResult.txt", coef2.dot(C2.T))
        np.savetxt("OriginalData.txt", 2.0*coef_fit[:, 0])
        #Res = coef2.dot(C2.T) - 2.0*coef_fit[:, 0]
        #np.savetxt("Residuals-M-Carbon.txt", Res)
        #print(Res)
    elif flag_se == "s":
        if coef3.shape[0] == np.linalg.matrix_rank(coef3):
            C3 = np.linalg.solve(coef3, 2.0*coef_fit[:, 1])
        else:
            C3 = np.linalg.lstsq(coef3, 2.0*coef_fit[:, 1], rcond = None)[0]
        #C2 = np.linalg.solve(coef2[[0, 1, 6], [0, 1, 6]], 2.0*coef_fit[[0, 1, 6], 0])
        C2 = np.linalg.lstsq(coef2, coef_fit[:, 0], rcond = None)[0]
    return (C2, C3)

def get_cijall(coef_fit, coef, Ord=3, flag_se="e"):
    '''
    Get the elastic constants of any order
    Parameters
    ----------
        coef_fit: np.ndarray
            The fitted coefficients of strain-energy/strain-stress equations
        coef: elastic3rd.symmetry.symmetry.CoefStr
            coef.coef2, coef.coef3, ...
        Ord: int
            The order of elastic constants
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
    Return
    ------
        C: elastic3rd.post.post.CCOEF
            The elastic constants, C.C2, C.C3, ...
    '''
    C = CCOEF(Ord)
    for i in range(2, int(Ord + 1)):
        a = 0
        if flag_se == "e":
            a = math.factorial(i)
        elif flag_se == "s":
            a = math.factorial(i - 1)
        NumCoef = eval('coef.coef' + str(i) + '.shape[0]')
        RankCoef = eval('np.linalg.matrix_rank(coef.coef' + str(i) + ')')
        if NumCoef == RankCoef:
            exec('C.C' + str(int(i)) + '= np.linalg.solve(coef.coef' + str(i) + ', a*coef_fit[:, i - 2])')
        else:
            exec('C.C' + str(int(i)) + '= np.linalg.lstsq(coef.coef' + str(i) + ', a*coef_fit[:, i - 2], rcond = None)[0]')
    return C

def get_cij_2nd(coef_2nd, coef2):
    '''
    Get the SOECs
    Parameters
    ----------
        coef_2nd: np.ndarray
            The fitting parameters of SOECs
        coef2: np.ndarray
            The analytical coefficient matrix of SOECs
    Return
    ------
        Cij2: np.ndarray
            SOECs
    '''
    Cij2 = np.linalg.lstsq(coef2, 2.0*coef_2nd, rcond = None)[0]
    return Cij2

def get_coef_2nd(s, e, V0):
    '''
    Get the fitting parameters of SOECs (for strain-energy method) (Not used)
    Parameters
    ----------
        s: np.ndarray
            strain
        e: np.ndarray
            energy
        V0: float
            volume
    Return
    ------
        coef_2nd: np.ndarray
            The fitted parameter of SOECs in strain-energy euqations
    '''
    eVpmol2GPa = 160.21719175
    (m, n) = e.shape
    coef_2nd = np.zeros((n, 1))
    for i in range(0, n):
        e0 = e[int((m-1)/2)][i]
        ei = e[:, i]
        ei = (ei - e0)/V0*eVpmol2GPa
        (coefi, pcovi) = esfit.esfit_2nd(s, ei)
        coef_2nd[i] = coefi
    return coef_2nd

def get_coef_nd(s, e, V0, flag_se, P0, acc=2):
    '''
    Get the coef by numerical differentiation method 
    (TODO: validate for the strain-energy method, probably not work)
    Parameters
    ----------
        s: np.ndarray
            strain
        e: np.ndarray
            stress
        V0: float
            volume
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
        P0: list/float
            The stress of un-deformed structure
        acc: int
            The accurace, 2, 4, 6 or 8
    Return
    ------
        coef_fit: np.ndarray
            The fitted coefficients of strain-stress equations
    '''
    if abs(s[0]) != abs(s[-1]):
        raise ValueError("For numerical differentiation, the strain must be symmetrical")
    flag_se = flag_se.lower()
    eVpmol2GPa = 160.21719175
    (m, n) = e.shape
    coef_fit = np.zeros((n, 2))
    for i in range(0, n):
        ei = e[:, i]
        e0 = esfit.get_e0(s, ei, flag_se)
        #e0 = e[n_d][i]  
        logger.info(e0)      
        logger.info(ei)
        if flag_se == "e":
            ei = (ei - e0)/V0*eVpmol2GPa
        elif flag_se == "s":
            ei = ei #(ei - e0)
        #TODO: give a general expression, ref: http://web.media.mit.edu/~crtaylor/calculator.html
        if acc == 2:
            #-1, 0, 1
            coef_fit[i, 0] = (ei[-1] - ei[0])/(2. * s[-1])
            coef_fit[i, 1] = (ei[-1] + ei[0] - 2.*P0[i])/(2. * s[-1] ** 2)
        elif acc == 4:
            #-2, -1, 0, 1, 2
            coef_fit[i, 0] = (8.*(ei[-2] - ei[1]) - (ei[-1] - ei[0]))/(12. * s[-2])
            coef_fit[i, 1] = (16.*(ei[-2] + ei[1]) - (ei[-1] + ei[0]) - 30.*P0[i])/(24. * s[-2] ** 2)
        elif acc == 6:
            #-3, -2, -1, 0, 1, 2, 3
            coef_fit[i, 0] = (45.*(ei[-3]-ei[2]) - 9.*(ei[-2]-ei[1]) + (ei[-1]-ei[0]))/(60. * s[-3])
            coef_fit[i, 1] = (270.*(ei[-3]+ei[2]) - 27.*(ei[-2]+ei[1]) + 2.*(ei[-1]+ei[0]) - 490.*P0[i])/(360.*s[-3]**2)
        elif acc == 8:
            #-4, -3, -2, -1, 0, 1, 2, 3, 4
            coef_fit[i, 0] = (672.*(ei[-4]-ei[3]) - 168.*(ei[-3]-ei[2]) + 32.*(ei[-2]-ei[1]) - 3.*(ei[-1]-ei[0]))/(840. * s[-4])
            coef_fit[i, 1] = (8064.*(ei[-4]+ei[3]) - 1008.*(ei[-3]+ei[2]) + 128.*(ei[-2]+ei[1]) - 9.*(ei[-1]+ei[0]) - 14350.*P0[i])/(10080.*s[-4]**2)
    return coef_fit

def get_coef(s, e, V0, flag_se, flag=3):
    '''
    Get the fitted coefficients by fitting
    Parameters
    ----------
        s: np.ndarray
            strain
        e: np.ndarray
            stress/energy
        V0: float
            volume
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
        flag: int
            The order of the polynomial
    Return
    ------
        coef_fit: np.ndarray
            The fitted coefficients of strain-stress equations

    '''
    #s: strain, e:energy
    flag_se = flag_se.lower()
    eVpmol2GPa = 160.21719175
    (m, n) = e.shape
    if flag > 3:
        coef_fit = np.zeros((n, flag - 1))
    else:
        coef_fit = np.zeros((n, 2))
    #n_d = int((m-1)/2)
    #s = np.delete(s, n_d)
    for i in range(0, n):
        ei = e[:, i]
        e0 = esfit.get_e0(s, ei, flag_se)
        #e0 = e[n_d][i]  
        logger.info(e0)      
        logger.info(ei)
        if flag_se == "e":
            ei = (ei - e0)/V0*eVpmol2GPa
        elif flag_se == "s":
            ei = ei #(ei - e0)
        if flag > 2:
            #print(len(s))
            #print(len(ei))
            (coefi, pcovi) = esfit.esfit(s, ei, flag_se, flag)
        else:
            s2 = s
            s2[n_d] = 1
            if flag == 1:
                e2 = ei/s2/s2
                s2 = np.delete(s2, n_d)
                e2 = np.delete(e2, n_d)
                (coefi, pcovi) = esfit.esfit(s2, e2, flag_se, flag)
            elif flag == 2:
                e2 = ei/s2
                #s2 = np.delete(s2, n_d)
                #e2 = np.delete(e2, n_d)
                s2[n_d] = 0
                (coefi, pcovi) = esfit.esfit(s2, e2, flag_se, flag)
        coef_fit[i, :] = coefi
    #if flag == 4:
    #    coef_fit = np.delete(coef_fit, -1, 1)
    return coef_fit

def read_e(EEnergy = "EEnergy.txt"):
    '''
    Read the energy or stress result
    Parameter
    ---------
        EEnergy: str
            The name of the energy/stress file
    Return
    ------
        e: np.ndarray
            The result of stress or energy
    '''
    e = np.loadtxt(EEnergy)
    return e


def post_mode(V0=1, Flag_Sep=1, Flag_Fig=1, Flag_Data=1, Flag_Ord=3, EEnergy="EEnergy.txt", INPUT="INPUT.e3rd", 
    STRAINMODE="STRAINMODE", Num_Str=1, FigName=None, DataFileName=None):
    '''
    Post process by specify strain modes
    Parameters
    ----------
        V0: float
            The volumn of the crystal, only required for Flag_SE = "e"
        Flag_Sep: 0 or 1
            0 for all as a whole, 1 for separate strains
        Flag_Fig: 0 or 1
            0 for don't show the fitting figures, 1 for show
        Flag_Data: 0 or 1
            1 for saving the higher order effect data, 0 for not
        Flg_Ord: int
            2-9, the order of polyfit used in fitting
        EEnergy: str
            The file contain the energy or stress, used in Elastic3rd
        INPUT: str
            The file contain the input parameters, used in Elastic3rd
        STRAINMODE: str
            The file contain the strain mode, used in Elastic3rd
        Num_Str: int
            The number of strain which was taken as a whole, only work for Flag_Sep = 1
        FigName: str
            The name of the figure, only work for Flag_Fig = 1
        DataFileName: str
            The name of the data file, only work for Flag_Data = 1
    Return
    ------
        C: elastic3rd.post.post.CCOEF
            The elastic constants, C.C2, C.C3, ...
    '''
    StrainIn = esutils.read_strainmode(STRAINMODE)
    (CrystalType, Ord, flag_se, StrainList) = get_post_param(INPUT)
    if not DataFileName:
        ParaIn = esutils.read_input(INPUT)
        DataFileName = ParaIn['BaseName']
    #print(StrainList)
    E = read_e(EEnergy)
    if Flag_Sep:
        m = StrainIn.shape
        m = m[0]
        k = m // Num_Str
        if m % Num_Str != 0:
            m = Num_Str * k
            logger.info('The number of strains can not be divided by Num_str. The extra strains are neglected.')
            StrainIn = StrainIn[0:m, :]
        for i in range(0, k):
            Ei = E[:, 6*Num_Str*i:6*Num_Str*(i+1)]
            Straini = np.zeros((Num_Str, 6))
            Straini[0:Num_Str, :] = StrainIn[Num_Str*i:Num_Str*(i+1)]
            C = post_single(StrainList/100., Ei, Straini, V0, Flag_Fig, Flag_Data, Flag_Ord, INPUT, FigName=FigName, DataFileName=DataFileName)
            logger.info(C.C2)
            logger.info(C.C3)
    else:
        C = post_single(StrainList/100., E, StrainIn, V0, Flag_Fig, Flag_Data, Flag_Ord, INPUT, FigName=FigName, DataFileName=DataFileName)
    return C

def post_single(x, E, StrainIn, V0, Flag_Fig=1, Flag_Data=1, Flag_Ord=3, INPUT="INPUT.e3rd", FigName=None, DataFileName=None):
    '''
    Post process for a single strain mode
    Parameters
    ----------
        x: np.ndarray
            strain
        E: np.ndarray
            Energy or Stress
        StrainIn: 2D list
            strain mode e.g. [[1, 2, -0.5, 1, 1, 0]]
        V0: float
            The volumn of the crystal, only required for Flag_SE = "e"
        Flag_Fig: 0 or 1
            0 for don't show the fitting figures, 1 for show
        Flag_Data: 0 or 1
            1 for saving the higher order effect data, 0 for not
        Flg_Ord: int
            2-9, the order of polyfit used in fitting
        INPUT: str
            The file contain the input parameters, used in Elastic3rd
        FigName: str
            The name of the figure, only work for Flag_Fig = 1
        DataFileName: str
            The name of the data file, only work for Flag_Data = 1
    Return
    ------
        C: elastic3rd.post.post.CCOEF
            The elastic constants, C.C2, C.C3, ...
    '''
    (CrystalType, Ord, flag_se, StrainList) = get_post_param(INPUT)
    C = CCOEF(Ord)
    if flag_se == "e":
        Cij_mode, coef_e, StrainMode = essym.CoefForSingleMode(CrystalType, Ord, StrainIn)
    elif flag_se == "s":
        Cij_mode, coef_e, StrainMode = symstr.CoefForSingleMode(CrystalType, Ord, StrainIn)
    coef_fit = get_coef(x, E, V0, flag_se, Flag_Ord)
    if Flag_Fig == 1:
        esfit.multiesplot(x, E, coef_fit, flag_se, Flag_Ord, V0, FigName=FigName)
    if Flag_Data == 1:
        esfit.get_data_effect_higherorder(x, E, coef_fit, flag_se, Flag_Ord, V0=V0, DataFileName=DataFileName)
    #print coef_fit
    #for i in range(2, Ord + 1):
    #    exec('coef' + str(i) + '= coef_e.coef' + str(i))
    #coef2 = coef_e.coef2
    #coef3 = coef_e.coef3
    #n is the column of coef3, which is equal to the number independent 3rd elastic constant in current crystal type
    n = eval('coef_e.coef' + str(int(Ord)) + '.shape[1]')
    SMRank = eval('np.linalg.matrix_rank(coef_e.coef' + str(int(Ord)) + ')')
    if SMRank < n:
        logger.info(SMRank)
        logger.info(n)
        logger.info('The strain mode is not full-rank.')
        #C2 = np.zeros((1, coef2.shape[1]))
        #C3 = np.zeros((1, n))
    else:
        C = get_cijall(coef_fit, coef_e, int(Ord), flag_se)
    return C

def post(V0, Flag_Fig = 1, Flag_Ord = 3, EEnergy = "EEnergy.txt", INPUT = "INPUT.e3rd", FigName=None):
    '''
    An old version of post function, not used now
    '''
    #StrainMode = read_strainmode(STRAINMODE)
    (CrystalType, Ord, flag_se, StrainList) = get_post_param(INPUT)
    E = read_e(EEnergy)
    if flag_se == "e":
        coef_e, StrainMode = essym.gen_strain_mode(CrystalType, Ord)
    elif flag_se == "s":
        coef_e, StrainMode = symstr.gen_strain_mode(CrystalType, Ord)
    coef_fit = get_coef(StrainList/100., E, V0, flag_se, Flag_Ord)
    #print coef_fit
    coef2 = coef_e.coef2
    coef3 = coef_e.coef3
    (C2, C3) = get_cij(coef_fit, coef2, coef3, flag_se)
    if Flag_Fig == 1:
        esfit.multiesplot(StrainList/100., E, coef_fit, flag_se, Flag_Ord, V0, FigName=FigName)
    return (C2, C3)

def get_post_param(INPUT="INPUT.e3rd"):
    '''
    Get post parameter from INPUT file
    Parameter
    ---------
        INPUT: str
            The name of INPUT file
    Return
    ------
        CrystalType: str
            The symmetry of the crystal
        Ord: int
            The order of elastic constant
        flag_se: str
            The strain-energy method(e) or strain-stress method(s)
        StrainList: list
            The strain list
    '''
    ParaIn = esutils.read_input(INPUT)
    flag_se = ParaIn['FlagSE'].lower()
    CrystalType = ParaIn['CrystalType']
    Ord = ParaIn['Ord']   
    StrainList = esutils.gen_strain_list(ParaIn)
    #print(StrainList)
    return (CrystalType, Ord, flag_se, StrainList)

class CCOEF:
    """This is the structure for coefficients. The attaches is coef + i, 
    where i is 2 to Ord
    Take Ord = 3 as an example, there are two attaches, coef2 and coef3"""
    def __init__(self, Ord = 3):
        for i in range(2, int(Ord) + 1):
            exec('self.C' + str(i) + ' = []')    
