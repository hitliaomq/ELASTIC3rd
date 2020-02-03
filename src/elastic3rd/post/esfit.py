#!python
#

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import math

def esfit_2nd(x, y):
    (coef, pcov) = curve_fit(esfun_2nd, x, y)
    return (coef, pcov)

def esfun_2nd(x, c2):
    y = c2 * x ** 2
    return y

def esfit_3rd(x, y):
    (coef, pcov) = curve_fit(esfun_3rd, x, y)
    return (coef, pcov)

def esfun_3rd(x, c2, c3):
    y = c2 * x ** 2 + c3 * x ** 3
    return y

## FITTING
def esfit(x, y, flag_se = "e", flag_ord = 3):
    flag_se = flag_se.lower()
    if flag_se == "e":
        coef, pcov = eval('curve_fit(esfun_energy_' + str(flag_ord) + ', x, y)')
        '''
        if flag_ord == 1:
            (coef, pcov) = curve_fit(esfun_energy_linear, x, y)
        elif flag_ord == 2:
            (coef, pcov) = curve_fit(esfun_energy_squa, x, y)
        elif flag_ord == 3:
            (coef, pcov) = curve_fit(esfun_energy_cubic, x, y)
        elif flag_ord == 4:
            (coef, pcov) = curve_fit(esfun_energy_quartic, x, y)
        elif flag_ord == 5:
            (coef, pcov) = curve_fit(esfun_energy_qui, x, y)
        '''
    elif flag_se == "s":
        pass
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
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 + c8 * x ** 8
    return y
def esfun_energy_9(x, c2, c3, c4, c5, c6, c7, c8, c9):
    y = c2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 + c5 * x ** 5 + c6 * x ** 6 + c7 * x ** 7 + c8 * x ** 8 + c9 * x ** 9
    return y

## PLOT
def yfitfun(x, c, flag_se = "e", flag_ord = 3):
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
        [c2, c3] = c
    flag_se = flag_se.lower()
    if flag_se == "e":
        #exec('')
        if flag_ord == 1:
            yfit = esfun_energy_1(x, c2, c3)
        elif flag_ord == 2:
            yfit = esfun_energy_2(x, c2, c3)
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

def esplot(x, y, coef, V0, flag_se = "e", flag_ord = 3):
    eVpmol2GPa = 160.21719175
    n = 100
    xmin = np.amin(x)
    xmax = np.amax(x)
    #xstep = (xmax - xmin)/n
    xfit = np.linspace(xmin, xmax, n)
    if flag_se == "e":
        yfit = V0/eVpmol2GPa*yfitfun(xfit, coef, flag_se, flag_ord)
    elif flag_se == "s":
        pass
    plotori = plt.plot(x, y, '*', label = 'original data')
    plotfit = plt.plot(xfit, yfit, 'r', label = 'fitting')
    plt.show()

def multiesplot(s, e, coef, flag_se, flag_ord, V0):
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
