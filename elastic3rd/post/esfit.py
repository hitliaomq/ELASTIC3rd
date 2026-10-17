#!python
#

from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from elastic3rd.esutils import logger

def esfit_2nd(x, y):
    '''
    Fitting y = a*x**2, and return the coefficients a and corresponding covariance
        Used for fitting the SOECs using strain-energy method
    '''
    (coef, pcov) = curve_fit(esfun_2nd, x, y)
    return (coef, pcov)

def esfun_2nd(x, c2):
    '''
    Define the function y = c2*x**2
    '''
    y = c2 * x ** 2
    return y

def esfit_3rd(x, y):
    '''
    Fitting y = a*x**2 + b*x**3, and return the coefficients a and corresponding covariance
        Used for fitting the TOECs using strain-energy method
    '''
    (coef, pcov) = curve_fit(esfun_3rd, x, y)
    return (coef, pcov)

def esfun_3rd(x, c2, c3):
    '''
    Define the function y = c2*x**2 + c3*x**3
    '''
    y = c2 * x ** 2 + c3 * x ** 3
    return y

## FITTING
def esfit(x, y, flag_se="e", flag_ord=3, sigma=None):
    '''
    Fitting the function for HOECs according to order and it will print the fitting result and correspondint errors
    Parameter
    ---------
        x, y: np.ndarray
            The strain(x) and energy/stress (y)
        flag_se: str
            The flag for the method, strain-energy method(e), strain-stress method(s)
        flag_ord: int
            The order of the elastic constant
        sigma: float
            Not used currently
    Return
    ------
        coef: np.ndarray
            The coefficients of the fitting, the number of elements is determined by the order
        pcov: np.ndarray
            The covariance of the fitting
    '''
    flag_se = flag_se.lower()
    #print(y)
    if flag_se == "e":
        coef, pcov = eval('curve_fit(esfun_energy_' + str(flag_ord) + ', x, y)')
        residuals = y - esfun_energy(x, coef)
        #residuals = eval('y - esfun_energy_' + str(flag_ord) + '(x, coef)')
    elif flag_se == "s":
        coef, pcov = eval('curve_fit(esfun_stress_' + str(flag_ord) + ', x, y)')
        residuals = y - esfun_stress(x, coef)
        #residuals = eval('y - esfun_stress_' + str(flag_ord) + '(x, coef)')
    else:
        raise Exception("Incorrect value for keyword (flag_SE)")
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    sigma_ab = np.sqrt(np.diagonal(pcov))
    logger.info("Fitting Result: ")
    logger.info(coef)
    logger.info("Error: ")
    logger.info(sigma_ab)
    logger.info("R^2 : " + str(r_squared))
    return (coef, pcov)

## FUNCTIONS
#for strain-energy method
#y:energy, x:strain
'''
def esfun(x, c, flag_ord):
    y = 0
    for i in range(0, flag_ord):
        exec('y = y + c[i] * x ** i')
    return y
'''
def esfun_energy(x, c):
    '''
    The general strain-energy function. The order is determined by the number of the elements in C (Ord=N+1)
    Parameters
    ----------
        x: float
            The strain
        c: list
            The elastic constants. [C2, C3, C4, ...]
    Return
    ------
        y: float
            The energy
    '''
    y = 0
    for i in range(0, len(c)):
        y = y + c[i] * (x ** (i + 2))
    return y

def esfun_energy_n(x, *args):
    '''
    The general strain-energy function. The order is determined by the number of inputs
    Parameters
    ----------
        x: float
            The strain
        args: several inputs
            The first is the C2, the second is C3, ...
    Return
    ------
        y: float
            The energy
    '''
    return esfun_energy(x, args)


#The following functions are the strain-energy function in different orders
def esfun_energy_3(x, c2, c3):
    y = c2 * x ** 2 + c3 * x ** 3
    return y
def esfun_energy_2(x, c2, c3):
    y = c2 * x + c3 * x ** 2
    return y
def esfun_energy_1(x, c2, c3):
    y = c2 + c3 * x
    return y
def esfun_energy_4(x, c2, c3, c4):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4
    return y
def esfun_energy_5(x, c2, c3, c4, c5):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5
    return y
def esfun_energy_6(x, c2, c3, c4, c5, c6):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6
    return y
def esfun_energy_7(x, c2, c3, c4, c5, c6, c7):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7
    return y
def esfun_energy_8(x, c2, c3, c4, c5, c6, c7, c8):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 \
        + c8 * x ** 8
    return y
def esfun_energy_9(x, c2, c3, c4, c5, c6, c7, c8, c9):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 \
        + c8 * x ** 8 + c9 * x ** 9
    return y
def esfun_energy_10(x, c2, c3, c4, c5, c6, c7, c8, c9, c10):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 \
        + c8 * x ** 8 + c9 * x ** 9 + c10 * x ** 10
    return y
def esfun_energy_11(x, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 \
        + c8 * x ** 8 + c9 * x ** 9 + c10 * x ** 10 + c11 * x ** 11
    return y

#for strain-stress method
#y:stress, x:strain

def esfun_stress(x, c):
    '''
    The general strain-stress function. The order is determined by the number of the elements in C (Ord=N+1)
    Parameters
    ----------
        x: float
            The strain
        c: list
            The elastic constants. [C2, C3, C4, ...]
    Return
    ------
        y: float
            The energy
    '''
    y = 0
    for i in range(0, len(c)):
        y = y + c[i] * (x ** (i + 1))
    return y

def esfun_stress_n(x, *args):
    '''
    The general strain-stress function. The order is determined by the number of input
    Parameters
    ----------
        x: float
            The strain
        args: several inputs
            The first is C2, The second one is C3, ...
    Return
    ------
        y: float
            The energy
    '''
    return esfun_stress(x, args)

#The following functions are the strain-stress function in different order
def esfun_stress_2(x, c2, c3):
    y = c2 + c3 * x
    return y
def esfun_stress_3(x, c2, c3):
    y = c2 * x + c3 * x ** 2
    return y

def esfun_stress_4(x, c2, c3, c4):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3
    return y
def esfun_stress_5(x, c2, c3, c4, c5):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4
    return y
def esfun_stress_6(x, c2, c3, c4, c5, c6):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5
    return y
def esfun_stress_7(x, c2, c3, c4, c5, c6, c7):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6
    return y
def esfun_stress_8(x, c2, c3, c4, c5, c6, c7, c8):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 \
        + c8 * x ** 7
    return y
def esfun_stress_9(x, c2, c3, c4, c5, c6, c7, c8, c9):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 \
        + c8 * x ** 7 + c9 * x ** 8
    return y
def esfun_stress_10(x, c2, c3, c4, c5, c6, c7, c8, c9, c10):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 \
        + c8 * x ** 7 + c9 * x ** 8 + c10 * x ** 9
    return y
def esfun_stress_11(x, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 \
        + c8 * x ** 7 + c9 * x ** 8 + c10 * x ** 9 + c11*x**10
    return y
def esfun_stress_12(x, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 \
        + c8 * x ** 7 + c9 * x ** 8 + c10 * x ** 9 + c11*x**10 + c12*x**11
    return y
def esfun_stress_13(x, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13):
    y = c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 \
        + c8 * x ** 7 + c9 * x ** 8 + c10 * x ** 9 + c11*x**10 + c12*x**11 + c13*x**12
    return y

## PLOT
def yfitfun(x, c, flag_se = "e", flag_ord = 3):
    '''
    General function for calculating energy/stress according to strain and elastic constants
    Parameters
    ----------
        x: float
            strain
        c: list
            the list of elastic constants
        flag_se: str
            strain-energy method (e) or strain-stress method (s)
        flag_ord: int
            The order of elastic constants
    Return
    ------
        yfit: float
            The energy/stress of current strain and elastic constants
    '''
    if flag_ord == 4:
        [c2, c3, c4] = c
    elif flag_ord == 5:
        [c2, c3, c4, c5] = c
    elif flag_ord == 6:
        [c2, c3, c4, c5, c6] = c
    elif flag_ord == 7:
        [c2, c3, c4, c5, c6, c7] = c
    elif flag_ord == 8:
        [c2, c3, c4, c5, c6, c7, c8] = c 
    elif flag_ord == 9:
        [c2, c3, c4, c5, c6, c7, c8, c9] = c 
    elif flag_ord == 10:
        [c2, c3, c4, c5, c6, c7, c8, c9, c10] = c
    elif flag_ord == 11:
        [c2, c3, c4, c5, c6, c7, c8, c9, c10, c11] = c
    elif flag_ord == 12:
        [c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12] = c
    elif flag_ord == 13:
        [c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13] = c
    elif flag_ord == 14:
        [c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14] = c
    elif flag_ord == 15:
        [c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15] = c
    else:
        [c2, c3] = c
    flag_se = flag_se.lower()
    input_str = []
    for i in range(0, len(c)):
        input_str.append("c" + str(2+i))
    input_str = ",".join(input_str)
    #print(input_str)
    if flag_se == "e":
        #exec('')
        if flag_ord == 1:
            yfit = esfun_energy_1(x, c2, c3)
        elif flag_ord == 2:
            yfit = esfun_energy_2(x, c2, c3)
        else:
            yfit = eval("esfun_energy_" + str(flag_ord) + "(x, " + input_str + ")")
    elif flag_se == "s":
        if flag_ord == 2:
            yfit = esfun_stress_2(x, c2, c3)
        else:
            yfit = eval("esfun_stress_" + str(flag_ord) + "(x, " + input_str + ")")
    return yfit

def esplot(x, y, coef, V0=None, flag_se = "e", flag_ord = 3):
    '''
    Plot the strain-energy or strain-stress points and its fitting result
    Parameters
    ----------
        x,y: np.ndarray
            The strain-energy or strain-stress points
        coef: list
            The fitting results
        V0: float
            The volume
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
        flag_ord: int
            The order of elastic constant
    Return
    ------
        This function will plot the figure and return the fitted value
    '''
    eVpmol2GPa = 160.21719175
    n = 201
    xmin = np.amin(x)
    xmax = np.amax(x)
    #xstep = (xmax - xmin)/n
    xfit = np.linspace(xmin, xmax, n)
    if flag_se == "e":
        yfit = V0/eVpmol2GPa*yfitfun(xfit, coef, flag_se, flag_ord)
    elif flag_se == "s":
        yfit = yfitfun(xfit, coef, flag_se, flag_ord)
    #data = np.zeros((n, 2))
    #data = np.concatenate((xfit.T, yfit.T), axis=0)
    #np.savetxt('data-HfCN-C0' + '.txt', yfit)
    plotori = plt.scatter(x, y)
    plotfit = plt.plot(xfit, yfit, label = 'fitting')
    #print(x)
    #print(y)
    #plt.savefig()
    #plt.show()
    return yfit

def effect_higher_order(x, y, coef, V0=None, flag_se = "e", flag_ord = 3):
    '''
    Get the first N effect
    Parameters
    ----------
        x,y: np.ndarray
            The strain-energy or strain-stress points
        coef: 1xN list
            The first N parameters in the fitting results
        V0: float
            The volume
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
        flag_ord: int
            The order of elastic constant
    Return
    ------
        data: np.ndarray
            The effect of first N effect
    '''
    eVpmol2GPa = 160.21719175
    n = 201
    xmin = np.amin(x)
    xmax = np.amax(x)
    #xstep = (xmax - xmin)/n
    xfit = np.linspace(xmin, xmax, n)
    data = xfit
    for i in range(len(coef)):
        coefi = coef[0 : i + 1]
        if i == 0:
            yfit = coefi * xfit
        else:
            yfit = yfitfun(xfit, coefi, flag_se, i + 2)
        if flag_se == "e":
            yfit = V0/eVpmol2GPa * yfit
        data = np.vstack((data, yfit))
    return data.T

def get_e0(s, e, flag_se):
    '''
    Get the energy/stress at initial structure
    Parameters
    ----------
        s: np.ndarray
            strain
        e: np.ndarray
            stress or energy
        flag_se: str
            strain-energy (e) or strain-stress (s) method
    Return
    ------
        e0: float
            The stress/energy of un-deformed structure
    '''
    n_d = np.where(np.abs(s) < 1e-5)
    #print(n_d)
    if n_d:
        #not empty
        n_d = n_d[0]
        e0 = e[n_d]
        #print(e0)
    else:
        #empty
        n_min = np.where(np.abs(s) == np.min(np.abs(s)))
        if len(n_min) > 1:
            n_min = np.where(s == np.min(s[n_min]))
        n_min = n_min[0]
        if n_min > 0:
            x1 = s[n_min]
            x2 = s[n_min - 1]
            y1 = e[n_min]
            y2 = e[n_min - 1]
            if flag_se == "e":
                #if is energy, the data around 0 is parabolic
                e0 = (x1 ** 2 * y2 - x2 ** 2 * y1)/(x1 ** 2 - x2 ** 2)
            elif flag_se == "s":
                #if is stress, the data around 0 is linear
                e0 = (x1 * y2 - x2 * y1)/(x1 - x2)
        else:
            e0 = np.min(s)
    return e0

def multiesplot(s, e, coef, flag_se, flag_ord=3, V0=None, FigName=None):
    '''
    Plot all s-e lines
    Parameters
    ----------
        s: 1D np.ndarray
            strain
        e: 2D np.ndarray
            energy/stress, each column is a line
        coef: 1xN list
            The first N parameters in the fitting results
        V0: float
            The volume
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
        flag_ord: int
            The order of elastic constant
        FigName: str
            The file name of the saved fig
                If the FigName is None (by default), show the figure
                Else, save the figure
    Return
    ------
        None
    '''
    if not FigName:
        matplotlib.use('Agg')
    (m, n) = e.shape
    n_d = int((m-1)/2)
    data = np.zeros((102, n))
    for i in range(0, n):
        coefi = coef[i, :]
        logger.info(coefi)
        ei = e[:, i]
        e0 = get_e0(s, ei, flag_se)
        if not e0:
            e0 = 0
        ei = ei - e0
        #e0 = e[n_d][i]
        if flag_ord > 2:
            s2 = s
            e2 = ei
        else:
            s2 = s
            s2[n_d] = 1
            if flag_ord == 1:
                e2 = ei/s2/s2
                s2 = np.delete(s2, n_d)
                e2 = np.delete(e2, n_d)
            elif flag_ord == 2:
                e2 = ei/s2
                #s2 = np.delete(s2, n_d)
                #e2 = np.delete(e2, n_d)
                s2[n_d] = 0
        yfit = esplot(s2, e2, coefi, V0, flag_se, flag_ord)
        #data[:, i] = yfit
        #np.savetxt('data-C-fit-Ord3-8' + '.txt', data)
    if FigName:
        plt.savefig(FigName, dpi=300)
        plt.cla()
    else:
        plt.show()

def get_data_effect_higherorder(s, e, coef, flag_se='s', flag_ord=3, V0=None, DataFileName=None):
    '''
    Get the each effect by considering n parameter
    Parameters
    ----------
        s, e: np.ndarray
            The strain-energy or strain-stress points
        coef: 1xN list
            The first N parameters in the fitting results
        V0: float
            The volume
        flag_se: str
            The flag for strain-stress method(s) or strain-energy method(e)
        flag_ord: int
            The order of elastic constant
        DataFileName: str
            The file name of the saved data, test by default
    Return
    ------
        None
    '''
    if not DataFileName:
        DataFileName = "test"
    (m, n) = e.shape
    for i in range(0, n):
        coefi = coef[i, :]
        logger.info(coefi)
        ei = e[:, i]
        e0 = get_e0(s, ei, flag_se)
        if not e0:
            e0 = 0
        ei = ei - e0        
        data = effect_higher_order(s, ei, coefi, V0=V0, flag_se=flag_se, flag_ord=flag_ord)
        np.savetxt(DataFileName + "-" + flag_se + "-" + str(flag_ord) + "-" + str(i) + '.txt', data)

'''
def mytestfun(x, a, b):
    y = 0
    #if flag == 1:
    y = a * x + b
    #elif flag == 2:
    #    y = a * x ** 2 + b * x
    return y

xdata = np.linspace(0, 4, 50)
y = mytestfun(xdata, 2.5, 1.3)
np.random.seed(1729)
y_noise = 0.2 * np.random.normal(size=xdata.size)
ydata = y + y_noise
plt.plot(xdata, ydata, 'b-', label='data')
plt.show()
(coef, pcov) = curve_fit(mytestfun, xdata, ydata)
print coef
'''

'''
#FOR TEST
x = np.arange(-0.05, 0.075, 0.025)
y = np.array([0.19287, 0.04706, 0, 0.04824, 0.17438])
(coef, pcov) = esfit(x, y)
coefall = np.zeros((3, 2))
print coef
coefall[0, :] = coef
print coefall
esplot(x, y, coef, 0)
#print x
#print y
'''
