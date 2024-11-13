import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from eigensystem import power_method

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

N_MAX = 50
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def test_power_method_1():
    '''
    Testing power method program by finding biggest eigenvalue of a stupidly simple matrix
    '''

    mat = np.diag([5, -1, 3])
    mat[-1][0] = 2
    print(mat)

    for i in np.geomspace(1, N_MAX, num = 20, dtype = int):
        plt.scatter(i, power_method(mat, i), marker = "x", s = 50, color = colors[2], label = "")

    print("Best approx for biggest eigenval:\t", power_method(mat, N_MAX))

    plt.hlines(5, 0, N_MAX, color = colors[0], label = "True value")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.ylabel("Iteration")
    plt.legend()
    plt.show()


def test_power_method_2():
    '''
    Testing power method program by finding biggest eigenvalue of a less stupidly simple matrix
    '''

    mat = np.array([[4, -1j, 2], [1j, 2, 2 + 7j], [2, 2 - 7j, -2]], dtype = complex)
    print(mat)

    # NOTE: Assuming eigenvalues are real:
    for i in np.geomspace(1, N_MAX, num = 20, dtype = int):
        plt.scatter(i, power_method(mat, i), marker = "x", s = 50, color = colors[2], label = "")

    print("Best approx for biggest eigenval:\t", power_method(mat, N_MAX))

    plt.hlines(8.45, 0, N_MAX, color = colors[0], label = "True value")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.ylabel("Iteration")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    test_power_method_1()
    test_power_method_2()