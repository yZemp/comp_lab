import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sympy

import sys
sys.path.append("..")
from matrices.matrix_utils import solve_linear_system
from interp_direct_monomial import interp_simple
from ex_1 import interp_newton
from ex_2 import chebyshev_nodes

from Polynomial_classes import Polynomial
from my_interplib import *

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

INTERVAL = (1, 4)
# INTERVAL = (-50, 50)
DELTA = 4
ORDER = 3
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]



# def _shift_row_right(mat, row, K):
#     n = mat.shape[1]  # Row length
#     if K > n : raise ValueError("Tried shifting a row for more than the length of the row")
#     mat[row] = np.concatenate((mat[row, - K :], mat[row, : - K]))


def _create_cubic_spline_mat(x, delta, order = 3, func = runge):
    '''
    Creates and returns the matrix for quadratic spline coefficients
    as well as the vector of ys.
    
    delta = number of sample points (including extrema)
    '''

    if order < 1: raise ValueError("The fuck bro")
    if order > 3: raise ValueError("Calm down")

    intervals = delta - 1
    n = order + 1 # Number of terms in every polynomial

    print("Delta:\t", delta)
    print("Order:\t", order)

    # Init mat and known_vector
    mat = np.zeros((n * intervals, n * intervals))
    known_vector = np.zeros(n * intervals)

    print("Fixing passage through points")

    for i in range(intervals): # Loop over intervals
        for k in range(n): # Loop over number of terms in polynomial
            print(i, k)
            
            # ALWAYS TWO ROWS (imposing to pass through two points)
            mat[i * intervals][k + i * n] = np.power(x[i], k)
            mat[i * intervals + 1][k + i * n] = np.power(x[i + 1], k)

            known_vector[i * intervals] = runge([x[i]])
            known_vector[i * intervals + 1] = runge([x[i + 1]])

            print(mat, "\n", known_vector)
            # exit()


    print("Fixing derivatives")

    offset = len(x)
    for i in range(1, intervals): # Loop over internal points
        for o in range(1, order): # Loop over orders of derivation
            xi = sympy.symbols('x')
            alpha = sympy.symbols('alpha')
            monomial = np.power(xi, alpha)
            deriv_o = sympy.diff(monomial, xi, o)
            
            for k in range(n): # Loop over monomials
                print(i, o, k)
                
                print(deriv_o)
                deriv_o_value = deriv_o.subs({alpha: k, xi: x[i]})
                
                mat[i + offset + o][k] = deriv_o_value
                mat[i + offset + o][k + n] = - deriv_o_value


            print(mat, "\n", known_vector)
            # exit()
            


    return mat, known_vector


def nth_spline(x, f, delta, order = 3, interval = INTERVAL):
    '''
    Perform piecewise interpolation of the runge function, given a specific sample
    with a polynomial of nth order

    We need to find a_i, b_i, c_i for every interval, so that
        f_i = a_i + b_i * x_i + c_i * x_i^2 + d_i * x_i^3 + ... + N_i * x_i^N
    '''

    if order < 1: raise ValueError("Order should be at least one.")

    mat, f_system = _create_nth_spline_mat(x, delta, order, func = runge)

    # Creating known array of f_i
    # FIXME: Generalize f_system to nth order 
    # Should be done, check if this works
    # f_system = np.array([[f[i], f[i + 1]] + [0] * (order - 1) for i in range(len(f) - 1)]).flatten()
    # print(len(f_system))
    print(mat, "\n", f_system)

    # Finding array of a_i, b_i, c_i
    interp = sp.linalg.solve(mat, f_system)
    # interp, _, _ = solve_linear_system(mat, f_system)
    print("Interp:\t", interp)


    print("Is the solution consistent:\t", np.allclose(np.dot(mat, interp), f_system))
    # print(np.dot(mat, interp))
    # print(f_system)

    # Creating polynomial for every interval
    # FIXME: Generalize f_system to nth order 
    # Should be done, check if this works
    polynomials = [Polynomial(interp[i:i + order + 1]) for i in range(0, len(interp), order + 1)]
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

    order = ORDER

    fig, ax = plt.subplots(1)
    # fig.suptitle("Uniform sampling")
    # plt.ylim(-1, 2)
    fig.set_size_inches(20, 10, forward = True)

    sample = np.linspace(*INTERVAL, num = DELTA)
    f = [runge(x) for x in sample]
    print(sample)
    arrx, arry, interp = nth_spline(sample, f, DELTA, order)
    # sample = np.linspace(*INTERVAL, num = 12)
    # sample, arrx, data, f = main(sample, 12, splines[0], color_index = 1)


    ax.scatter(sample, f, color = (.1, .1, .1), marker = "x", label = "Samples")
    ax.plot(arrx, arry, color = colors[0], alpha = 1, label = f"Spline {len(sample)} order: {order}")
    ax.plot(arrx, runge(arrx), c = (.3, .3, .3), alpha = .1, linewidth = 10, label = "Runge function")
    ax.legend()


    plt.savefig(f"interp_graphs/nth_spline_runge")
    plt.show()