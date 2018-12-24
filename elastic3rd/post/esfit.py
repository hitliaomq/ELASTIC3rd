#!python
#

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import math

#e:energy, s:strain
def esfit(x, y):
    (coef, pcov) = curve_fit(esfun, x, y)
    return (coef, pcov)

def esfun(x, c2, c3):
    y = c2 * x ** 2 + c3 * x ** 3
    return y

def essolve():
    pass

def esplot(x, y, coef, saveflag):
    n = 100
    xmin = np.amin(x)
    xmax = np.amax(x)
    xstep = (xmax - xmin)/n
    c2 = coef[0]
    c3 = coef[1]
    xfit = np.arange(xmin, xmax + xstep, xstep)
    yfit = esfun(xfit, c2, c3)
    plotori = plt.plot(x, y, '*', label = 'original data')
    plotfit = plt.plot(xfit, yfit, 'r', label = 'fitting')
    plt.show()

#'''
#FOR TEST
x = np.arange(-0.05, 0.075, 0.025)
y = np.array([0.19287, 0.04706, 0, 0.04824, 0.17438])
(coef, pcov) = esfit(x, y)
coefall = np.zeros((3, 2))
print coef
coefall[0, :] = coef
print coefall
#esplot(x, y, coef, 0)


#print x
#print y
#'''