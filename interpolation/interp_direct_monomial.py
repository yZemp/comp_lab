import numpy as np
import matplotlib.pyplot as plt
import random as rand
import scipy.stats as sts

import sys
sys.path.append("../")
from matrices.matrix_utils import get_inverse, mat_vec_prod
from interpolation.Polynomial_classes import Polynomial

# Direct approach using Vandermonde (monomial) matrix

np.set_printoptions(suppress = True, precision = 3)

def get_vandermonde_mat(arr):
    '''
    Given an array, this function returns its respective vandermonde matrix
    '''
    return np.array([[np.power(x, i) for i in range(len(arr))] for x in arr])

def interp_simple(x, f):
    vandermonde_mat = get_vandermonde_mat(x)
    anti_vandermonde = get_inverse(vandermonde_mat)
    a = np.dot(anti_vandermonde, f)
    print(vandermonde_mat, "\n")
    print(np.linalg.inv(vandermonde_mat), "\n")
    print(anti_vandermonde, "\n")

    lnsp = np.linspace(-5, 35, 1_000)

    # Creation of polynomial to be evaluated
    poly = Polynomial(a)
    print(a, "\n", poly)
    # Evaluating polynomial in all of lnsp
    y = [poly.evaluate(el) for el in lnsp]
    return lnsp, y


if __name__ == "__main__":
    # Starting data
    glob_x = np.array([0, 10, 15, 20, 22.5, 30]) # i = 0, 1, 2, 3, 4, 5
    glob_f = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.67])

    # Splitting data:
    start = 2
    end = 4
    x = glob_x[start:end]
    f = glob_f[start:end]

    plt.scatter(glob_x, glob_f, marker = "x", color = (.9, .1, .1), label = "Fixed points")

    lnsp, y = interp_simple(x, f)
    plt.plot(lnsp, y, color = (.2, .3, .9), label = "Simplest fit")


    plt.legend()
    plt.show()