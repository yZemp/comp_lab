import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from matrices.matrix_utils import solve_linear_system
from interpolation.ex_2 import chebyshev_nodes

from interpolation.Polynomial_classes import Polynomial
from interpolation.my_interplib import *

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

INTERVAL = (-1, 1)
# INTERVAL = (-50, 50)
DELTA = 11
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def _create_quadratic_spline_mat(x, delta):
    '''
    Creates and returns the matrix for quadratic spline coefficients (it's a quasi-block matrix)
    
    delta = number of sample points (including extrema)
    '''

    print("Delta:\t", delta)

    # Building a diagonal block matrix
    blocks = [[[np.power(x[k + i], j) if k + i < len(x) else 0. for j in range(3)] for i in range(3)] for k in range(delta - 1)]
    # print(blocks)
    mat = sp.linalg.block_diag(*blocks)
    
    # Fixing derivatives rows
    for r in range(2, len(mat), 3):
        for i in range(r - 2, r - 2 + 6):
            # Subbing only if in bounds
            if i < len(mat[r]):
                sign = 1 if i < (r + 1) else - 1
                # print(sign, i, r)
                # print(mat)
                if i % 3 == 0: mat[r][i] = 0
                if i % 3 == 1: mat[r][i] = sign
                if i % 3 == 2: mat[r][i] = 2 * sign * x[(r + 1) // 3]

    return np.array(mat)


def quadratic_spline(x, f, delta, interval = INTERVAL):
    '''
    Perform piecewise interpolation of the runge function, given a specific sample
    with a polynomial of second order

    We need to find a_i, b_i, c_i for every interval, so that
        f_i = a_i + b_i * x_i + c_i * x_i^2
    '''

    mat = _create_quadratic_spline_mat(x, delta)

    # Creating known array of f_i
    f_system = np.array([[f[i], f[i + 1], 0] for i in range(len(f) - 1)]).flatten()
    # print(len(f_system))
    print(mat, "\n", f_system)

    # Finding array of a_i, b_i, c_i
    # interp = sp.linalg.solve(mat, f_system)
    interp, _, _ = solve_linear_system(mat, f_system)
    print("Interp:\t", interp)


    # print("Is the solution consistent:\t", np.allclose(np.dot(mat, interp), f_system))
    # print(np.dot(mat, interp))
    # print(f_system)

    # Creating polynomial for every interval
    polynomials = [Polynomial(interp[i:i + 2 + 1]) for i in range(0, len(interp), 2 + 1)]
    # print(len(interp), len(polynomials))
    # print(polynomials)

    arrx = np.linspace(*interval, num = 1_000)
    arry = [polynomials[max(0, next(i for i, xbar in enumerate(x) if xbar >= el) - 1)].evaluate(el) for el in arrx]

    return arrx, arry, interp



if __name__ == "__main__":
    '''
    Perform a piecewise interpolation.

    NOTE: delta is always == number_of_intervals + 1
    '''

    fig, ax = plt.subplots(2)
    # fig.suptitle("Uniform sampling")
    # plt.ylim(-1, 2)
    fig.set_size_inches(20, 10, forward = True)

    sample = np.linspace(*INTERVAL, num = DELTA)
    f = [runge(x) for x in sample]
    print(sample)
    arrx, arry, interp = quadratic_spline(sample, f, DELTA)
    # sample = np.linspace(*INTERVAL, num = 12)
    # sample, arrx, data, f = main(sample, 12, splines[0], color_index = 1)


    ax[0].scatter(sample, f, color = (.1, .1, .1), marker = "x", label = "Samples (dispari)")
    ax[0].plot(arrx, arry, color = colors[2], alpha = 1, label = f"Spline {len(sample)}")
    ax[0].plot(arrx, runge(arrx), c = (.3, .3, .3), alpha = .1, linewidth = 10, label = "Runge function")
    ax[0].legend()


    sample = np.linspace(*INTERVAL, num = DELTA + 1)
    cheby = chebyshev_nodes(DELTA)
    f = [runge(x) for x in sample]
    print(sample)
    arrx, arry, interp = quadratic_spline(sample, f, DELTA + 1)
    # sample = np.linspace(*INTERVAL, num = 12)
    # sample, arrx, data, f = main(sample, 12, splines[0], color_index = 1)


    ax[1].scatter(sample, f, color = (.1, .1, .1), marker = "x", label = "Samples (pari)")
    ax[1].plot(arrx, arry, color = colors[2], alpha = 1, label = f"Spline {len(sample)}")
    ax[1].plot(arrx, runge(arrx), c = (.3, .3, .3), alpha = .1, linewidth = 10, label = "Runge function")
    ax[1].legend()


    plt.savefig(f"interp_graphs/quadratic_spline_runge")
    plt.show()