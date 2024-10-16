import numpy as np
import matplotlib.pyplot as plt
import random as rand
import scipy.stats as sts

import sys
sys.path.append("~/Documents/Programming/comp_lab/")
# from matrices.matrix_utils import get_inverse
from interpolation.Polynomial_classes import *

# Direct approach using Newton polynomials


def _interp(arr):
    '''
    Should receive array of array of f_i
    Returns array of divided differences
    '''
    if len(arr) == 1: return arr

    # FIXME

    # Here x is global data
    temp = np.array([(arr[-1][i + 1] - arr[-1][i]) / (x[i + 1] - x[i]) for i in range(0, (arr[-1]) - 1)], dtype = np.float128)
    turbomat = np.append(arr, _interp(temp))
    return turbomat


# Init data
x = np.array([0, 10, 15, 20, 22.5, 30]) # i=0,1,2,3,4,5
f = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.67])


if __name__ == "__main__":

    # Derivation of interpolation coefficients

    # Init starting as array of array of f_i
    # This will become (in the next line) the turbomatrix
    starting = np.array([f])
    # Derive complete turbomatrix
    turbomatrix = _interp(starting)
    
    print(turbomatrix)

    # plt.legend()
    # plt.show()
    
    
    # plt.plot(x, y, c = (.6, 1 - .2 * i, .8), label = f"Order: {i}")