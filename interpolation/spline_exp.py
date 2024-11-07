import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from matrices.matrix_utils import solve_linear_system
from interp_direct_monomial import interp_simple
from ex_1 import interp_newton

from Polynomial_classes import Polynomial
from my_interplib import *

#######################################################################
# VARS

INTERVAL = (-1, 1)
DELTA = 11
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def spline(x, f):
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

    return interp


def main(sample, label, delta):
    '''
    Perform piecewise interpolation of the runge function, given a specific sample
    '''

    f = [runge(x) for x in sample]
    interp = spline(sample, f)
    # Creating polynomial for every interval
    polynomials = [Polynomial(interp[i:i + 2]) for i in range(0, len(interp) - 2, 2)]
    # print(len(polynomials))

    arrx = np.linspace(*INTERVAL, num = 1_000)
    data = [polynomials[min(max(next(i for i, xbar in enumerate(sample) if xbar >= x) - 1, 0), delta - 3)].evaluate(x) for x in arrx]
    # data = [max(next(i for i, xbar in enumerate(sample) if xbar >= x) - 1, 0) for x in arrx]
    # print(data)
    # exit()

    return sample, arrx, data, f



if __name__ == "__main__":
    '''
    Perform a piecewise interpolation.
    '''

    # TODO: Unsure how to generate the interval
    # intervals = np.linspace(*INTERVAL, num = int((INTERVAL[1] - INTERVAL[0]) / DELTA))

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
        uniform = np.linspace(*INTERVAL, num = n)
        sample, arrx, data, f = main(uniform, "Uniform", n)
    
        if n == start: ax.plot(arrx, runge(arrx), c = (.3, .3, .3), alpha = .1, linewidth = 10, label = "Runge function")
        ax.plot(arrx, data, color = (.1 * n / 5, 1 - .2 * n / 10, .8), alpha = 1, label = f"Spline {len(sample)}")

            

    plt.legend()
    plt.savefig(f"interp_graphs/spline")
    plt.show()