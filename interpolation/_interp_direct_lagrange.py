import numpy as np
import matplotlib.pyplot as plt
import random as rand
import scipy.stats as sts

import sys
sys.path.append("../")
from matrices.matrix_utils import get_inverse
from interpolation.Polynomial_classes import *

# Direct approach using Lagrange polynomials 

def get_interp_mat(arr):
    '''
    Given an array, this function returns its respective interpolation matrix
    In theis case it's the identity
    '''
    return np.array([[1. if i == j else 0. for j in range(len(arr))] for i in range(len(arr))], dtype = float)


if __name__ == "__main__":
    # Starting data
    x = np.array([0, 10, 15, 20, 22.5, 30]) # i=0,1,2,3,4,5
    f = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.67])

    # x = np.sort([100 * rand.random() for i in range(50)])
    # f = [el ** 3 - 2 * el + sts.norm.rvs(0, 1) for el in x]
    
    # plt.yscale("log")

    plt.scatter(x, f, marker = "x", c = "#040404", label = "Data")

    interp_mat = get_interp_mat(x)
    anti_mat = get_inverse(interp_mat)
    a = np.dot(anti_mat, f)


    lnsp = np.linspace(-5, 35, 1_000)

    # Creation of lagrangian polynomial
    poly = Lag_polynomial(a)
    # Evaluating polynomial in all of lnsp
    y = [poly.evaluate(el) for el in lnsp]
    plt.plot(lnsp, y, c = "#a30404", label = "Fit")


    plt.legend()
    plt.show()
    
    
    # plt.plot(x, y, c = (.6, 1 - .2 * i, .8), label = f"Order: {i}")