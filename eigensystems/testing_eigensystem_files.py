import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from eigensystem import power_method, power_method_all, inverse_power_method

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

N_MAX = 50
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def separator(n):
    print("\n")
    print("".join(["#" for _ in range(100)]), f"\n# {n}")


def plot_img(z):
    plt.hlines(np.real(z), 0, N_MAX, color = colors[0], label = "True value (Re)")
    plt.hlines(np.imag(z), 0, N_MAX, color = colors[1], label = "True value (Im)")

def plot_img_V2(z):
    plt.hlines(np.abs(z), 0, N_MAX, color = colors[0], label = "True value mod")



def test_power_method_1():
    '''
    Testing power method program by finding biggest eigenvalue of a stupidly simple matrix
    '''

    separator(1)

    mat = np.diag([5, -1, 3])
    mat[-1][0] = 2
    print(mat)

    for i in np.geomspace(1, N_MAX, num = 20, dtype = int):
        plt.scatter(i, power_method(mat, i), marker = "x", s = 50, color = colors[2], label = "")

    print("Best approx for biggest eigenval:\t", power_method(mat, N_MAX))

    big_eigenval = np.sort(np.linalg.eigvals(mat))[-1]
    plot_img(big_eigenval)
    print("True biggest eigenvalue:\t", big_eigenval)
    # plt.xscale("log")
    # plt.yscale("log")
    plt.ylabel("Iteration")
    plt.legend()
    plt.show()


def test_power_method_2():
    '''
    Testing power method program by finding biggest eigenvalue of a less stupidly simple matrix
    '''

    separator(2)

    mat = np.array([[4, -1j, 2], [1j, 2, 2 + 7j], [2, 2 - 7j, -2]], dtype = complex)
    print(mat)

    for i in np.geomspace(1, N_MAX, num = 20, dtype = int):
        plt.scatter(i, power_method(mat, i), marker = "x", s = 50, color = colors[2], label = "")

    print("Best approx for biggest eigenval:\t", power_method(mat, N_MAX))

    big_eigenval = np.sort(np.linalg.eigvals(mat))[-1]
    plot_img(big_eigenval)
    print("True biggest eigenvalue:\t", big_eigenval)
    # plt.xscale("log")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.ylabel("Iteration")
    plt.legend()
    plt.show()


def test_power_method_3():
    '''
    Testing power method program by finding biggest eigenvalue of a not so simple matrix
    '''

    separator(3)

    mat = np.array([[5, -1, 3, 0], [-8, 26, 1903, 666], [2, 0, 7, -26], [0, 4, 4, -44]])
    print(mat)

    for i in np.geomspace(1, N_MAX, num = 20, dtype = int):
        plt.scatter(i, power_method(mat, i), marker = "x", s = 50, color = colors[2], label = "")

    print("Best approx for biggest eigenval:\t", power_method(mat, N_MAX))

    # NOTE: Biggest eigenvalue in abs
    big_eigenval = np.sort(np.linalg.eigvals(mat))[0]
    plot_img(big_eigenval)
    print("True biggest eigenvalue:\t", big_eigenval)
    # plt.xscale("log")
    # plt.yscale("log")
    plt.ylabel("Iteration")
    plt.legend()
    plt.show()


def test_power_method_all_1():
    '''
    Using the power method program to find all eigenvalues of a stupidly simple matrix
    '''

    separator(4)

    np.random.seed(0)
    mat = np.dot(np.random.randn(3, 3).T, np.random.randn(3, 3))
    print(mat)

    eigenvalues, eigenvectors = power_method_all(mat, N_MAX)

    print(np.sort(eigenvalues))
    print(np.sort(np.linalg.eigvals(mat)))


def test_power_method_all_2():
    '''
    Using the power method program to find all eigenvalues of a less stupidly simple matrix
    '''

    separator(5)

    mat = np.array([[4, -1j, 2], [1j, 2, 2 + 7j], [2, 2 - 7j, -2]], dtype = complex)
    print(mat)

    eigenvalues, eigenvectors = power_method_all(mat, N_MAX)

    print(np.sort(eigenvalues))
    print(np.sort(np.linalg.eigvals(mat)))



def test_inverse_power_method_1():
    '''
    Testing inverse power method program by finding smallest eigenvalue of a stupidly simple matrix
    '''

    separator(6)

    mat = np.array([[4, -1j, 2], [1j, 2, 2 + 7j], [2, 2 - 7j, -2]], dtype = complex)
    print(mat)

    for i in np.geomspace(1, 10, num = 20, dtype = int):
        plt.scatter(i, inverse_power_method(mat, i, shifted = False), marker = "x", s = 50, color = colors[2], label = "")

    print("Best approx for smallest eigenval:\t", inverse_power_method(mat, 3, shifted = False))

    # NOTE: convergence is stupidly fast
    small_eigenval = np.sort([abs(z) for z in np.linalg.eigvals(mat)])[0]
    plot_img_V2(small_eigenval)
    print("True smallest eigenvalue:\t", small_eigenval)
    # plt.xscale("log")
    # plt.yscale("log")
    plt.ylabel("Iteration")
    plt.legend()
    plt.show()


def converg_study():
    '''
    Testing inverse power method program by finding smallest eigenvalue of a stupidly simple matrix
    '''

    separator(7)

    mat = np.array([[4, -1j, 2], [1j, 2, 2 + 7j], [2, 2 - 7j, -2]], dtype = complex)
    print(mat)

    N_MAX = 100

    big = []
    small = []
    small_shifted = []

    arri = []
    big_T = power_method(mat, N_MAX)
    small_T = inverse_power_method(mat, N_MAX, shifted = False)
    small_shifted_T = inverse_power_method(mat, N_MAX, shifted = True)

    for i in np.geomspace(1, N_MAX, num = 50, dtype = int):
        arri.append(i)
        big.append(np.abs(power_method(mat, i) / big_T - 1))
        small.append(np.abs(inverse_power_method(mat, i, shifted = False) / small_T - 1))
        small_shifted.append(np.abs(inverse_power_method(mat, i, shifted = True) / small_shifted_T - 1))

    plt.figure(figsize = (10, 6))

    plt.plot(arri, big, color = colors[0], label = "Convergence power method")
    plt.plot(arri, small, color = colors[2], label = "Convergence inverse power method")
    plt.plot(arri, small_shifted, color = colors[1], linewidth = 5, alpha = .5, zorder = 0, label = "Convergence inverse power method shifted")

    plt.xlabel("Max iteration")
    plt.ylabel("Error")
    plt.yscale("log")

    plt.legend()
    plt.show()

if __name__ == "__main__":
    # test_power_method_1()
    # test_power_method_2()
    # test_power_method_3()
    # test_power_method_all_1()
    test_power_method_all_2()
    # test_inverse_power_method_1()
    # converg_study()