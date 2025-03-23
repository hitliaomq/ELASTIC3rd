#!python
#

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import math

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
def esfit(x, y, flag_se = "e", flag_ord = 3, ec_order=3):
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
    Return
    ------
        coef: np.ndarray
            The coefficients of the fitting, the number of elements is determined by the order
        pcov: np.ndarray
            The covariance of the fitting
    '''
    flag_se = flag_se.lower()
    if ec_order == 3:
        if flag_ord == 2:
            pass
    if flag_se == "e":
        if ec_order == 3:
            if flag_ord == 2:
                coef, pcov = curve_fit(esfun_energy_3_2, x, y)
            elif flag_ord == 1:
                coef, pcov = curve_fit(esfun_energy_3_1, x, y)
            else:
                coef, pcov = eval('curve_fit(esfun_energy_' + str(flag_ord) + ', x, y)')
        else:
            if flag_ord < ec_order:
                raise ValueError
            else:
                coef, pcov = eval('curve_fit(esfun_energy_' + str(flag_ord) + ', x, y)')
    elif flag_se == "s":
        pass
    return (coef, pcov)


#The following functions are the strain-energy equations by considering different higher order effect
def esfun_energy_3(x, c2, c3):
    y = c2 * x ** 2 + c3 * x ** 3
    return y
def esfun_energy_2(x, c2):
    y = c2 * x ** 2
    return y
def esfun_energy_3_2(x, c2, c3):
    y = c2 * x + c3 * x ** 2
    return y
def esfun_energy_3_1(x, c2, c3):
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
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 + c8 * x ** 8
    return y
def esfun_energy_9(x, c2, c3, c4, c5, c6, c7, c8, c9):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 + c8 * x ** 8 + c9 * x ** 9
    return y

## PLOT
def yfitfun(x, c, flag_se = "e", flag_ord = 3, ec_order=3):
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
    else:
        if ec_order == 3:
            [c2, c3] = c
        else:
            if flag_ord == 3:
                [c2, c3] = c
            elif flag_ord == 2:
                c2 = c
            else:
                raise ValueError
    flag_se = flag_se.lower()
    if flag_se == "e":
        #exec('')
        if flag_ord == 1:
            if ec_order == 3:
                yfit = esfun_energy_3_1(x, c2, c3)
            else:
                raise ValueError
        elif flag_ord == 2:
            if ec_order == 3:
                yfit = esfun_energy_3_2(x, c2, c3)
            else:
                yfit = esfun_energy_2(x, c2)
        elif flag_ord == 3:
            yfit = esfun_energy_3(x, c2, c3)
        elif flag_ord == 4:
            yfit = esfun_energy_4(x, c2, c3, c4)
        elif flag_ord == 5:
            yfit = esfun_energy_5(x, c2, c3, c4, c5)
        elif flag_ord == 6:
            yfit = esfun_energy_6(x, c2, c3, c4, c5, c6)
        elif flag_ord == 7:
            yfit = esfun_energy_7(x, c2, c3, c4, c5, c6, c7)
        elif flag_ord == 8:
            yfit = esfun_energy_6(x, c2, c3, c4, c5, c6, c7, c8)
        elif flag_ord == 9:
            yfit = esfun_energy_7(x, c2, c3, c4, c5, c6, c7, c8, c9)
    elif flag_se == "s":
        pass
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
    n = 100
    xmin = np.amin(x)
    xmax = np.amax(x)
    #xstep = (xmax - xmin)/n
    xfit = np.linspace(xmin, xmax, n)
    if flag_se == "e":
        yfit = V0/eVpmol2GPa*yfitfun(xfit, coef, flag_se, flag_ord, ec_order=flag_ord)
    elif flag_se == "s":
        pass
    plotori = plt.plot(x, y, '*', label = 'original data')
    plotfit = plt.plot(xfit, yfit, 'r', label = 'fitting')
    plt.show()

def multiesplot(s, e, coef, flag_se='e', flag_ord=3, V0=None):
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
    (m, n) = e.shape
    n_d = int((m-1)/2)
    for i in range(0, n):
        coefi = coef[i, :]
        print(coefi)
        ei = e[:, i] - e[n_d][i]
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
        esplot(s2, e2, coefi, V0, flag_se, flag_ord)


