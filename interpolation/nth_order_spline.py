import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sympy

import sys
sys.path.append("..")

from .Polynomial_classes import Polynomial
from .my_interplib import *

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

INTERVAL = (-20, 20)
# INTERVAL = (-50, 50)
DELTA = 30
ORDER = 3
LINE_WIDTH = 8
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

def FUNC(x):
    return np.cos(x) * x * sp.stats.norm.pdf(x, 0, 5)
    # return np.cos(x) * x


# def _shift_row_right(mat, row, K):
#     n = mat.shape[1]  # Row length
#     if K > n : raise ValueError("Tried shifting a row for more than the length of the row")
#     mat[row] = np.concatenate((mat[row, - K :], mat[row, : - K]))


def _create_cubic_spline_mat(x, f, delta, order = 3):
    '''
    Creates and returns the matrix for quadratic spline coefficients
    as well as the vector of ys.
    
    delta = number of sample points (including extrema)
    '''

    intervals = delta - 1
    n = order + 1 # Number of terms in every polynomial

    print("Delta:\t", delta)
    print("Order:\t", order)

    # Init mat and known_vector
    mat = np.zeros((n * intervals, n * intervals), dtype = complex)
    known_vector = np.zeros((n * intervals), dtype = complex)

    print("Fixing passage through points")

    for i in range(intervals): # Loop over intervals
        for k in range(n): # Loop over number of terms in polynomial
            # print(i, k)
            
            # ALWAYS TWO ROWS (imposing to pass through two points)
            mat[i * 2][k + i * n] = np.power(x[i], k)
            mat[i * 2 + 1][k + i * n] = np.power(x[i + 1], k)

            known_vector[i * 2] = f[i]
            known_vector[i * 2 + 1] = f[i + 1]

            # print(mat, "\n", known_vector)
            # exit()


    print("Fixing derivatives")

    offset = intervals * 2
    xi = sympy.symbols('x')
    alpha = sympy.symbols('alpha')
    for i in range(1, intervals): # Loop over internal points
        for o in range(1, order): # Loop over orders of derivation
            monomial = np.power(xi, alpha)
            deriv_o = sympy.diff(monomial, xi, o)
            
            for k in range(0, n): # Loop over monomials
                # print(i, o, k)
                
                # print(deriv_o)
                deriv_o_value = deriv_o.subs({alpha: k, xi: x[i]})
                
                mat[offset + (i - 1) * (order - 1) + o - 1][(i - 1) * n + k] = deriv_o_value
                mat[offset + (i - 1) * (order - 1) + o - 1][(i - 1) * n + k + n] = - deriv_o_value


            # print(mat, "\n", known_vector)
            # exit()

    print("Fixing remaining dof")
    
    if order >= 2:

        xi = sympy.symbols('x')
        alpha = sympy.symbols('alpha')
        for o in range(1, order): # Loop over orders of derivation
            monomial = np.power(xi, alpha)
            deriv_o = sympy.diff(monomial, xi, o)
            
            for k in range(0, n): # Loop over monomials
                deriv_o_value = deriv_o.subs({alpha: k, xi: x[-1]})
                mat[- order + 1][- n + k] = deriv_o_value

    if order >= 3:

        xi = sympy.symbols('x')
        alpha = sympy.symbols('alpha')
        for o in range(1, order): # Loop over orders of derivation
            monomial = np.power(xi, alpha)
            deriv_o = sympy.diff(monomial, xi, o)
            
            for k in range(0, n): # Loop over monomials
                deriv_o_value = deriv_o.subs({alpha: k, xi: x[0]})
                mat[- order + 2][k] = deriv_o_value


    # print("FINAL FUCKING MATRIX (and vector):\n")
    # print(mat, "\n", known_vector)


    return mat, known_vector


def nth_spline(x, f, delta, order = 3, interval = INTERVAL):
    '''
    Perform piecewise interpolation of the runge function, given a specific sample
    with a polynomial of nth order

    We need to find a_i, b_i, c_i for every interval, so that
        f_i = a_i + b_i * x_i + c_i * x_i^2 + d_i * x_i^3 + ... + N_i * x_i^N
    '''

    if order < 1: raise ValueError("Order should be at least one.")
    if order > 3: raise ValueError("Chill the fuck out bro.")

    mat, f_system = _create_cubic_spline_mat(x, f, delta, order = order)

    # Finding array of a_i, b_i, c_i
    interp = sp.linalg.solve(mat, f_system)
    # interp, _, _ = solve_linear_system(mat, f_system)
    # print("Interp:\t", interp)


    print("Is the solution consistent:\t", np.allclose(np.dot(mat, interp), f_system))
    # print(np.dot(mat, interp))
    # print(f_system)

    # Creating polynomial for every interval
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
    # myfunc = runge
    myfunc = FUNC

    fig, ax = plt.subplots(1)
    # fig.suptitle("Uniform sampling")
    # plt.ylim(-1, 2)
    fig.set_size_inches(20, 10, forward = True)

    sample = np.linspace(*INTERVAL, num = DELTA)
    f = [myfunc(x) for x in sample]
    print(sample)
    arrx, arry, interp = nth_spline(sample, f, DELTA, order)
    # sample = np.linspace(*INTERVAL, num = 12)
    # sample, arrx, data, f = main(sample, 12, splines[0], color_index = 1)


    ax.scatter(sample, f, color = (.1, .1, .1), marker = "x", label = "Samples")
    ax.plot(arrx, arry, color = colors[0], alpha = 1, label = f"Spline {len(sample)} order: {order}")
    ax.plot(arrx, myfunc(arrx), c = (.3, .3, .3), alpha = .1, linewidth = LINE_WIDTH, label = "True function")
    ax.legend()


    plt.savefig(f"interp_graphs/nth_spline_runge")
    plt.show()