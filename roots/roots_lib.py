import numpy as np
import matplotlib.pyplot as plt

MAX_ITER = 1e6
PREC = 1e-9

def bisect(a, b, f, prec = PREC, max_iter = MAX_ITER, _i = 0):
    '''
    Given a function and an interval returns the root inside said interval, assuming:
        the root exists
        the root is unique in said interval
        the function changes sign before and after the root
    '''
    
    if _i > max_iter or abs(a - b) <= prec: return a + (b - a) / 2
    else: _i += 1

    c = a + (b - a) / 2
    if f(a) * f(c) < 0: return bisect(a, c, f, prec, max_iter, _i)
    else: return bisect(c, b, f, prec, max_iter, _i)


def newton_raphson(x0, f, g, prec = PREC, max_iter = MAX_ITER, _i = 0, fix_oscillations = True):
    '''
    Given a function, its derivative and one initial guess returns the root, assuming:
        the root exists

    NOTE: A good initial guess is had after some bisection iterations
    '''

    # Exit conditions
    if _i > max_iter or abs(f(x0)) < prec: return x0 
    else: _i += 1

    x1 = x0 - f(x0) / g(x0)

    # Check and fix oscillation
    if fix_oscillations:
        if x1 - f(x1) / g(x1) == x0: x1 = (x1 + x0) / 2
        
    return newton_raphson(x1, f, g, prec, max_iter, _i, fix_oscillations = fix_oscillations)


def secant(x0, x1, f, prec = PREC, max_iter = MAX_ITER, _i = 0):
    '''
    Given a function and two initial guesses returns the root, assuming:
        the root exists
        the function is derivable

    NOTE: A good initial guess is had after some bisection iterations
    This function relies on the backward finite difference to approximate the first derivative
    '''

    # Exit conditions
    if _i > max_iter or abs(f(x1)) < prec: return x1 
    else: _i += 1

    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))

    return secant(x1, x2, f, prec, max_iter, _i)




    







if __name__ == "__main__":
    def f(x): return .5 + x - np.power(x, 2)

    root = bisect(.8, 1.6, f)

