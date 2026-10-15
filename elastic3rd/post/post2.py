#!python
#
#New post script for elastic3rd, this script will replace post.py in the future
#Features:
#   1. using polyfit in numpy
#   2. weighted fitting (all using numerical differentiation?)

import numpy as np
import math
from scipy.optimize import curve_fit
import elastic3rd.tools.nd_coef as esnd
import elastic3rd.symmetry.symmetry as essym
import elastic3rd.symmetry.symmetry_stress as symstr
import elastic3rd.esutils as esutils
from elastic3rd.esutils import logger

def esfun_energy(x, c):
    '''
    Function for strain-energy method with list input

    Parameter
    ---------
        x: list/np.ndarray
            strain
        c: list/np.ndarray
            The coefficients of strain-energy equation
    Return
    ------
        y: float
            The energy. 
    '''
    y = 0
    for i, ci in enumerate(c):
        y = y + ci * (x ** (i + 2)) / float(math.factorial(i + 2))
    return y

def esfun_energy_n(x, *args):
    '''
    Function for strain-energy method with multi inputs

    Parameter
    ---------
        x: list/np.ndarray
            strain
        *args: multi args
            The coefficients of the strain-energy equation
    Return
    ------
        Energy
    '''
    return esfun_energy(x, args)

def esfun_energy_low1(x, c):
    '''
    This function used to hand the strain-energy method, and the y(energy density) is divided by x(strain)
    for the higher-order effect is not obvious in energy but obvious in stress
    The input parameter is the same with esfun_energy
    '''
    y = 0
    for i, ci in enumerate(c):
        y = y + ci * (x ** (i + 1)) / float(math.factorial(i + 2))
    return y

def esfun_energy_low1_n(x, *args):
    '''
    Thw multi input parameters version of esfun_energy_low1
    '''
    return esfun_energy_low1(x, args)

def esfun_energy_low2(x, c):
    '''
    This function used to hand the strain-energy method, and the y(energy density) is divided by x^2(strain)
    for the higher-order effect is not obvious in energy but obvious in stress
    The input parameter is the same with esfun_energy
    '''
    y = 0
    for i, ci in enumerate(c):
        y = y + ci * (x ** i) / float(math.factorial(i + 2))
    return y

def esfun_energy_low2_n(x, *args):
    '''
    Thw multi input parameters version of esfun_energy_low2
    '''
    return esfun_energy_low2(x, args)

def esfun_stress(x, c):
    '''
    Function for strain-stress method with list input

    Parameter
    ---------
        x: list/np.ndarray
            strain
        c: list/np.ndarray
            The coefficients of strain-stress equations
    Return
    ------
        y: float
            Stress
    '''
    y = 0
    for i, ci in enumerate(c):
        y = y + ci * (x ** (i + 1)) / float(math.factorial(i + 1))
    return y

def esfun_stress_n(x, *args):
    '''
    The multi input parameters version of esfun_stress
    '''
    return esfun_stress(x, args)

def esfun_stress_low1(x, c):
    '''
    This function used to hand the strain-stress method, and the y(stress) is divided by x(strain)
    The input parameter is the same with esfun_stress
    '''
    y = 0
    for i, ci in enumerate(c):
        y = y + ci * (x ** i) / float(math.factorial(i + 1))
    return y

def esfun_stress_low1_n(x, *args):
    '''
    The multi input parameters version of esfun_stress_low1
    '''
    return esfun_stress_low1(x, args)

def yfitfun(x, c, flag_se="s", flag_ord=3, flag_low=0):
    '''
    General function for calculating energy/stress according to strain and elastic constants
    Parameters
    ----------
        x: np.ndarray
            Strain
        c: np.ndarray
            The coefficients of the strain-energy or strain-stress equations
        flag_se: str
            strain-energy method(e) or strain-stress method(s)
        flag_ord: int
            The order of the elastic constants
        flag_low: int
            Use the low-order method or not. for flag_se='s', it can be 0 or 1; for flag_se='e', it can be 0 or 1 or 2
    Return
    ------
        yfit: float
            The value of the strain-energy or strain-stress functions
    '''
    if flag_low == 0:
        if flag_se == 's':
            yfit = esfun_stress(x, c)
        elif flag_se == 'e':
            yfit = esfun_energy(x, c)
    elif flag_low == 1:
        if flag_se == 's':
            yfit = esfun_stress_low1(x, c)
        elif flag_se == 'e':
            yfit = esfun_energy_low1(x, c)
    elif flag_low == 2:
        if flag_se == 'e':
            yfit = esfun_energy_low2(x, c)
        else:
            raise ValueError('The flag_low for strain-stress method can be 0 or 1, current value is {}'.format(flag_low))
    else:
        raise ValueError('The value for flag_low should no more than 2. Current value is flag_low={}'.format(flag_low))
    return yfit

def get_fitting_report(yfit, ydata, pcov):
    '''
    Print the fitting report, including the fitting error and MSE and R2
    Parameters
    ----------
        yfit: list/np.ndarray
            The fitted value(energy or stress) of strain-(energy/stress) function
        ydata: list/np.ndarray
            The original data of energy or stress
        pcov: np.ndarray
            The covariance of the fitting
    Return
    ------
        perr: np.ndarray
            The error of the fitting
        R2: float
            The coefficient of determination
        MSE: float
            The Mean Squared Error
        SSE: float
            The sum of squares for error
        RMSE: float
            The root of MSE
    '''
    yfit = np.array(yfit)
    ydata = np.array(ydata)
    n = ydata.size
    # SSE（Sum of Squares for Error） = RSS(residual sum of squares)
    SSE = ((ydata - yfit)**2).sum()
    # MSE (Mean Squared Error) or called Standard Deviation
    MSE = SSE/n
    # RMSE is the root of MSE
    RMSE = np.sqrt(MSE)
    
    u = ydata.mean()
    # SST(Sum of Squares fortotal) = TSS(total sum of squares)
    SST = ((ydata - u)**2).sum()
    # SSR(Sum of Squares forregression) = ESS (explained sum of squares)
    SSR = SST - SSE
    #for coefficient of determination
    R2 = SSR/SST

    perr = np.sqrt(np.diagonal(pcov))

    logger.info('The error for the parameters is {}'.format(perr.tolist()))
    logger.info('The Mean Squared Error(MSE) is MSE={}'.format(MSE))
    logger.info('The coefficient of determination is R^2={}'.format(R2))
    return (perr, R2, MSE, SSE, RMSE)

def esfit(x, ys, V0=1, flag_se="s", flag_ord=3, weight=False, flag_low=0):
    '''
    Multi fittings with weight(weight=True) or without (weight=False)
    Parameters
    ----------
        x: np.ndarray
            strain
        ys: 2D np.ndarray
            energy or stress
        V0: float
            volume
        flag_se: str
            strain-energy method(e) or strain-stress method(s)
        flag_ord: int
            The order of elastic constant
        weight: bool
            Fitting with weight or not
        flag_low: int
            Using lower-order or not
    Return
    ------
        coef_fit: np.ndarray
            The coefficient of the strain-energy/stress equation by fitting method
        perr: np.ndarray
            The fitting error
        R2: np.ndarray
            The coefficients of determination
        MSE: np.ndarray
            The Mean Square Error
    '''
    eVpmol2GPa = 160.21719175
    (m, n) = ys.shape
    coef_fit, perr, R2, MSE = [], [], [], []
    for i in range(n):
        y = ys[:, i]
        if flag_se == 'e':
            y = y/V0*eVpmol2GPa
        if weight:
            (coefi, perri, R2i, MSEi) = esfit_weighted(x, y, flag_se=flag_se, flag_ord=flag_ord, flag_low=flag_low)
        else:
            (coefi, perri, R2i, MSEi) = esfit_noweight(x, y, flag_se=flag_se, flag_ord=flag_ord, flag_low=flag_low)
        coef_fit.append(coefi)
        perr.append(perri)
        R2.append(R2i)
        MSE.append(MSEi)
    return (np.array(coef_fit), np.array(perr), np.array(R2), np.array(MSE))

def esfit_noweight(x, y, flag_se="s", flag_ord=3, sigma=None, flag_low=0):
    '''
    Single fitting the equations without weight
    Parameters
    ----------
        x: np.ndarray
            strain
        y: np.ndarray
            energy density or stress
        flag_se: str
            strain-energy method(e) or strain-stress method(s)
        flag_ord: int
            The order of elastic constant
        sigma: float
            sigma in curve_fit. This parameter is only used when calling current function by esfit_weighted
        flag_low: int
            Using lower-order or not
    Return
    ------
        coef: np.ndarray
            The fitting result
        perr: np.ndarray
            The fitting error
        R2: np.ndarray
            The coefficients of determination
        MSE: np.ndarray
            The Mean Square Error
    '''
    flag_se = flag_se.lower()
    x0 = [1 for i in range(flag_ord-1)]
    if flag_low == 0:
        if flag_se == "e":
            coef, pcov = curve_fit(esfun_energy_n, x, y, p0=x0, sigma=sigma)
        elif flag_se == "s":
            coef, pcov = curve_fit(esfun_stress_n, x, y, p0=x0, sigma=sigma)
        else:
            raise ValueError("Incorrect value for keyword (flag_se)")
    elif flag_low == 1:
        y = [y[i]/xi if xi != 0 else 0 for i, xi in enumerate(x)]
        if flag_se == "e":
            coef, pcov = curve_fit(esfun_energy_low1_n, x, y, p0=x0, sigma=sigma)
        elif flag_se == "s":
            coef, pcov = curve_fit(esfun_stress_low1_n, x, y, p0=x0, sigma=sigma)
        else:
            raise ValueError("Incorrect value for keyword (flag_se)")
    elif flag_low == 2:
        y = [y[i]/(xi ** 2) if xi != 0 else 0 for i, xi in enumerate(x)]
        if flag_se == "e":
            coef, pcov = curve_fit(esfun_energy_low2_n, x, y, p0=x0, sigma=sigma)
        else:
            raise ValueError("Incorrect value for flag_low in strain-stress method (should be 0 or 1, currently is {})".format(flag_low))
    else:
        raise ValueError('The value for flag_low should no more than 2. Current value is flag_low={}'.format(flag_low))
    yfit = yfitfun(x, coef, flag_se=flag_se, flag_ord=flag_ord, flag_low=flag_low)
    logger.info('The fitting result is {}'.format(coef))
    (perr, R2, MSE, SSE, RMSE) = get_fitting_report(yfit, y, pcov)
    return (coef, perr, R2, MSE)

def esfit_weighted(x, y, flag_se='s', flag_ord=3, flag_low=0):
    '''
    Fitting the equations using weighted data points
        Note: The weight is taken from the formula of numerical differentiation in current x
              For different degree, the weight is different
    Parameters
    ----------
        x: np.ndarray
            strain
        y: np.ndarray
            energy density or stress
        flag_se: str
            strain-energy method(e) or strain-stress method(s)
        flag_ord: int
            The order of elastic constant
        flag_low: int
            Using lower-order or not
    Return
    ------
        coef: np.ndarray
            The fitting result
        perr: np.ndarray
            The fitting error
        R2: np.ndarray
            The coefficients of determination
        MSE: np.ndarray
            The Mean Square Error
    '''
    flag_se = flag_se.lower()
    coef, R2, MSE, perr = [], [], [], []
    for i in range(0, flag_ord-1):
        if flag_se == 'e':
            deg = i + 2
            w = esnd.gen_nd_coef_e(x, deg=deg)
        elif flag_se == 's':
            deg = i + 1
            w = esnd.gen_nd_coef(x, deg=deg)
        (x0, y0, sigma) = get_sigma_from_w(x, y, w[:, 0], deg=deg)
        logger.info(x0)
        logger.info(w[:, 0])
        (coefi, perri, R2i, MSEi) = esfit_noweight(x0, y0, flag_se=flag_se, flag_ord=flag_ord, sigma=sigma, flag_low=flag_low)
        coef.append(coefi[i])
        perr.append(perri[i])
        R2.append(R2i)
        MSE.append(MSEi)
    return (coef, perr, R2, MSE)

def get_sigma_from_w(x, y, w, deg=2):
    '''
    Get the sigma from weight, and if wi=0, remove corresponding xi and yi
    Parameters
    ----------
        x: np.ndarray
            strain
        y: np.ndarray
            stress or energy
        w: np.ndarray
            weight
        deg: int
            The order of the small strain
    Return
    ------
        x: np.ndarray
            strain
        y: np.ndarray
            stress or energy
        sigma: np.ndarray
            sigma in curve_fit. Here it equals to normalized 1/w
    '''
    tor = 1e-5
    h = min(abs(i) for i in x if abs(i) > tor)
    w = np.absolute(w)
    x = [x[i] for i, wi in enumerate(w) if wi > tor/(h**deg)]
    y = [y[i] for i, wi in enumerate(w) if wi > tor/(h**deg)]
    w = [wi for i, wi in enumerate(w) if wi > tor/(h**deg)]
    sigma = 1./np.array(w)
    sigma_max = np.max(sigma)
    sigma = [sigmai/sigma_max for sigmai in sigma]
    return (np.array(x), np.array(y), np.array(sigma))

def get_coef(x, ys, V0=1, flag_se='s', flag_ord=3):
    '''
    Get the coefficient by numerical differentiation
        Note: ys is mxn np.ndarray, m is the datapoints for each fitting, and n is the number of fittings
    Parameters
    ----------
        x: np.ndarray
            strain
        ys: 2D np.ndarray
            energy or stress
        V0: float
            volume
        flag_se: str
            strain-energy method(e) or strain-stress method(s)
        flag_ord: int
            The order of elastic constant
    Return
    ------
        coef_nd: np.ndarray
            The coefficient of the strain-energy/stress equation by numerical differentiation result
    '''
    eVpmol2GPa = 160.21719175
    (m, n) = ys.shape
    coef_nd = []
    for i in range(n):
        y = ys[:, i]
        if flag_se == 'e':
            y = y/V0*eVpmol2GPa
        coef_nd.append(get_coef_single(x, y, flag_se=flag_se, flag_ord=flag_ord))
    return np.array(coef_nd)

def get_coef_single(x, y, flag_se='s', flag_ord=3):
    '''
    Get the coef of nd's result
    Here x and y are 1d np.ndarray
    Parameters
    ----------
        x: np.ndarray
            strain
        ys: 2D np.ndarray
            energy or stress
        flag_se: str
            strain-energy method(e) or strain-stress method(s)
        flag_ord: int
            The order of elastic constant
    Return
    ------
        coef_nd: np.ndarray
            The coefficient of the strain-energy/stress equation by numerical differentiation result
    '''
    coef_nd = []
    if flag_se == 's':
        deg_start = 1
        deg_end = 3
    elif flag_se == 'e':
        deg_start = 2
        deg_end = 4
    for deg in range(deg_start, deg_end):
        if flag_se == 's':
            coef_nd.append(np.dot(y, esnd.gen_nd_coef(x, deg=deg))[0])
        elif flag_se == 'e':
            coef_nd.append(np.dot(y, esnd.gen_nd_coef_e(x, deg=deg))[0])
    return coef_nd

def get_cij(coef_fit, coef, flag_se='s', Ord=3):
    '''
    Get the elastic constant by the fitting parameters
    Parameters
    ----------
        coef_fit: np.ndarray
            The fitting parameters
        coef: elastic3rd.symmetry.symmetry.CoefStr
            coef.coef2, coef.coef3, ...
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
        Ord: int
            The order of elastic constants
    Return
    ------
        C: elastic3rd.post.post.CCOEF
            The elastic constants, C.C2, C.C3, ...
    '''
    flag_se = flag_se.lower()

    C = []
    for i in range(Ord-1):
        if flag_se.startswith('s'):
            pre_factor = math.factorial(Ord-i-1)
        else:
            pre_factor = math.factorial(Ord-i)
        coefi = eval('coef.coef{}'.format(Ord-i))
        if coefi.shape[0] == np.linalg.matrix_rank(coefi):
            Ci = np.linalg.solve(coefi, pre_factor * coef_fit[:, Ord-2-i]).tolist()
        else:
            Ci = np.linalg.lstsq(coefi, pre_factor * coef_fit[:, Ord-2-i], rcond = None)[0].tolist()
        C.append(Ci)
    return C

def get_post_param(INPUT = "INPUT.e3rd"):
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

def post_single(x, E, StrainIn, CoefMethod='nd', V0=1, Flag_Ord=3, INPUT="INPUT.e3rd", flag_low=0):
    '''
    Post for a single strain mode
    Parameters
    ----------
        x: np.ndarray
            strain
        E: np.ndarray
            Energy or Stress
        StrainIn: 2D list
            strain mode e.g. [[1, 2, -0.5, 1, 1, 0]]
        CoefMethod: str
            The method to get the coef
                nd: numerical differentiation
                fit: fitting without weighted
                fitw: fitting with weight
        V0: float
            The volumn of the crystal, only required for Flag_SE = "e"
        Flg_Ord: int
            2-9, the order of polyfit used in fitting
        INPUT: str
            The file contain the input parameters, used in Elastic3rd
        flag_low: int
            Using the lower-order equations (divide the strain)
    Return
    ------
        C: elastic3rd.post.post.CCOEF
            The elastic constants, C.C2, C.C3, ...
        Cerr: elastic3rd.post.post.CCOEF
            The error on elastic constants, Cerr.C2, Cerr.C3, ... (Note: only work for fitting method)
    
    '''
    CoefMethod = CoefMethod.lower()
    (CrystalType, Ord, flag_se, StrainList) = get_post_param(INPUT)
    if flag_se == "e":
        Cij_mode, coef_e, StrainMode = essym.CoefForSingleMode(CrystalType, Ord, StrainIn)
    elif flag_se == "s":
        Cij_mode, coef_e, StrainMode = symstr.CoefForSingleMode(CrystalType, Ord, StrainIn)
    if CoefMethod == 'nd':
        coef_fit = get_coef(x, E, V0=V0, flag_se=flag_se, flag_ord=Ord)
        #perr = np.zeros(len(coef_fit))
    elif CoefMethod == 'fit':
        (coef_fit, perr, R2, MSE) = esfit(x, E, V0=V0, flag_se=flag_se, flag_ord=Flag_Ord, weight=False, flag_low=flag_low)
    elif CoefMethod == 'fitw':
        (coef_fit, perr, R2, MSE) = esfit(x, E, V0=V0, flag_se=flag_se, flag_ord=Flag_Ord, weight=True, flag_low=flag_low)
    C = get_cij(coef_fit, coef_e, flag_se=flag_se, Ord=Ord)
    if CoefMethod == 'nd':
        return C
    else:
        Cerr = get_cij(perr, coef_e, flag_se=flag_se, Ord=Ord)
        return (C, Cerr)

