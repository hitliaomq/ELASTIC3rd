#!python
#
#Generate the coefficients of numerical differentiation
#
#Ref.http://web.media.mit.edu/~crtaylor/calculator.html
#    https://finite-difference-coefficients.nl/
#
#Part of Elastic3rd
#
# Mingqing Liao
# liaomq1900127@163.com

import numpy as np

def gen_int_x(x):
    h = min(abs(i) for i in x if abs(i) > 1e-5)
    return (h, [xi/h for xi in x])

def gen_coef_matrix(x, flagse='s'):
    A = []
    n = len(x)
    if flagse == 'e':
        n = len(x) + 1
    for i in range(n):
        A.append([xi**i for xi in x])
    return A

def gen_b(n_x, deg=2, flagse='s'):
    if flagse == 'e':
        n_x = n_x + 1
    deg = int(deg)
    b = np.zeros((n_x))
    b[deg] = np.math.factorial(deg)
    return b

def gen_nd_coef(x, deg=2):
    '''
    Parameters
    ----------
        x: list
            [0, 1, 2] means f(x), f(x+h), f(x+2h)
        deg: int
            The order of deriver
    '''
    n_x = len(x)
    (h, x) = gen_int_x(x)
    A = np.array(gen_coef_matrix(x))
    b = np.array(gen_b(n_x, deg=deg)).reshape((-1, 1))
    return np.dot(np.linalg.inv(A), b) / h**deg

def gen_nd_coef_e(x, deg=2):
    '''
    Get the numerical differentiation coefficients of energy method,
    Compared with stress method, the first-order term always equals to 0
    '''
    n_x = len(x)
    (h, x) = gen_int_x(x)
    A = np.array(gen_coef_matrix(x, flagse='e'))
    A = np.delete(A, 1, axis=0)
    b = np.array(gen_b(n_x, deg=deg, flagse='e')).reshape((-1, 1))
    b = np.delete(b, 1, axis=0)
    return np.dot(np.linalg.inv(A), b) / h**deg

