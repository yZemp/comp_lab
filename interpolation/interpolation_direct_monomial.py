import numpy as np
import matplotlib.pyplot as plt
import random as rand
import scipy.stats as sts

import sys
sys.path.append("../")
from matrices.matrix_utils import get_inverse, mat_vec_prod

# Direct approach using Vandermonde (monomial) matrix


np.set_printoptions(suppress=True,precision=3)


def get_vandermonde_mat(arr):
    '''
    Given an array, this function returns its respective vandermonde matrix
    '''
    return np.array([[np.power(x, i) for i in range(len(arr))] for x in arr], dtype = np.float128)
    

def polyn(param):
    '''
    Returns a polynomial of n-th order to be evaluated
    
    if param is an array: 
        param is array of coefficients
    
    if param is int: 
        polynomial order is param, and all coefficients are assumed to be 1

    '''

    if type(param) == int:
        if param == 0: return lambda x: 1

        def func(x):
            return np.sum([np.power(x, i) for i in range(0, param)], dtype = np.float128)
        
        return func

    if type(param) == np.ndarray:

        def func(x):
            return np.sum([el * np.power(x, i) for i, el in enumerate(param)])

        return func
    
    
    return -1


if __name__ == "__main__":
    # Starting data
    x = np.array([0, 10, 15, 20, 22.5, 30]) # i=0,1,2,3,4,5
    f = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.67])

    # x = np.sort([100 * rand.random() for i in range(50)])
    # f = [el ** 3 - 2 * el + sts.norm.rvs(0, 1) for el in x]
    
    # plt.yscale("log")

    plt.scatter(x, f, marker = "x", c = "#040404", label = "Data")

    vandermonde_mat = get_vandermonde_mat(x)
    anti_vandermonde = get_inverse(vandermonde_mat)
    a = np.dot(anti_vandermonde, f)


    lnsp = np.linspace(-5, 35, 1_000)

    # Creation of polynomial to be evaluated
    poly = polyn(a)
    # Evaluating polynomial in all of lnsp
    y = [poly(el) for el in lnsp]
    plt.plot(lnsp, y, c = "#a30404", label = "Fit")


    plt.legend()
    plt.show()
    
    
    # plt.plot(x, y, c = (.6, 1 - .2 * i, .8), label = f"Order: {i}")