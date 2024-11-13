import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from matrices.matrix_utils import solve_linear_system

from Polynomial_classes import Polynomial
from my_interplib import *

#######################################################################
# VARS

INTERVAL = (-1, 1)
DELTA = 11
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def linear_spline(x, f, delta = DELTA, interval = INTERVAL):
    '''
    Perform a piecewise interpolation given points with coordinate x, f

    We need to find a_i, b_i for every interval, so that
        f_i = a_i + b_i * x_i
    '''

    # Creating block matrix
    blocks = [[[1, x[i]], [1, x[i + 1]]] for i in range(len(x) - 1)]
    mat = sp.linalg.block_diag(*blocks)
    # print(mat)

    # Creating known array of f_i
    f_system = np.array([[f[i], f[i]] for i in range(len(f))]).flatten()[1:-1]
    # print(f_system)

    # Finding array of a_i, b_i
    interp, _, _ = solve_linear_system(mat, f_system)
    # print(interp)

    # Creating polynomial for every interval
    polynomials = [Polynomial(interp[i:i + 2]) for i in range(0, len(interp), 2)]
    # print(len(interp), len(polynomials))

    arrx = np.linspace(*interval, num = 1_000)
    arry = [polynomials[max(0, next(i for i, xbar in enumerate(x) if xbar >= el) - 1)].evaluate(el) for el in arrx]
    # print(data)
    # exit()

    return arrx, arry, interp


if __name__ == "__main__":
    '''
    Perform a piecewise interpolation.
    '''
    
    # Uniform sample

    fig, ax = plt.subplots(1)
    # fig.suptitle("Uniform sampling")
    plt.ylim(-1, 2)
    fig.set_size_inches(20, 10, forward = True)

    # ax.scatter(sample, f, color = (.1, .1, .1), label = "Data sample", s = 70, marker = "x")

    start = 6
    end = 50


    # Showing the convergence
    for i in np.geomspace(start, end, num = 10):
        n = int(i)
        # print(n)
        sample = np.linspace(*INTERVAL, num = n)
        f = [runge(x) for x in sample]
        arrx, data, interp = linear_spline(sample, f, n)
    
        if n == start: ax.plot(arrx, runge(arrx), c = (.3, .3, .3), alpha = .1, linewidth = 10, label = "Runge function")
        ax.plot(arrx, data, color = (.1 * n / 5, 1 - .2 * n / 10, .8), alpha = 1, label = f"Spline {len(sample)}")

            

    plt.legend()
    plt.savefig(f"interp_graphs/spline")
    plt.show()