import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from matrices.matrix_utils import solve_linear_system

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

N_MAX = 30
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def _rayleigh_quotient(A, yn):
    '''
    Returns the rayleigh quotient of a square matrix and a vector
    '''

    num = np.vdot(np.transpose(yn), np.dot(A, yn))
    den = np.vdot(yn, yn)
    return num / den

def power_method(mat, N = N_MAX, return_eigenvector = False):
    '''
    Finds and returns the biggest eigenvalue of a square matrix
    '''

    np.random.seed(0)
    # x0 = np.array([row[0] for row in mat])
    x0 = np.random.random(len(mat)) + 1j * np.random.random(len(mat))
    y_n = [x0]

    for _ in range(N):
        new_y = np.dot(mat, y_n[-1])
        new_y_normalized = new_y / np.linalg.norm(new_y)
        y_n.append(new_y_normalized)
    
    approx = _rayleigh_quotient(mat, y_n[-1])
    if return_eigenvector: return approx, y_n[-1]
    else: return approx


def _power_method_all(mat, N = N_MAX):
    '''
    Finds eigensystem of a matrix using deflation
    '''

    eigenvalues = []
    eigenvectors = []
    current_mat = np.copy(mat)

    print("Calculating eigenvalues...")

    for _ in range(len(mat)):

        eigenval, eigenvec = power_method(current_mat, N, return_eigenvector = True)
        eigenvec = eigenvec / np.linalg.norm(eigenvec)
        eigenvalues.append(eigenval)
        eigenvectors.append(eigenvec)
        
        current_mat = current_mat - eigenval * np.outer(eigenvec, eigenvec.conj())

    return eigenvalues, eigenvectors

def power_method_all(mat, N = N_MAX):
    '''
    Wrapper for _power_method_all
    Finds eigensystem of a matrix using deflation, shifting mat by a scalar.
    This way, we can make every eigenvalue positive without changing eigenvectors.
    Then we can shift the spectrum back
    '''

    first_eigenval = power_method(mat, N)
    shifted_mat = mat + 2 * first_eigenval * np.identity(len(mat))
    # print("Biggest eigenvalue:\t", first_eigenval)
    # print(shifted_mat)
    shifted_eigenvalues, eigenvectors = _power_method_all(shifted_mat, N)
    eigenvalues = shifted_eigenvalues - 2 * first_eigenval

    return eigenvalues, eigenvectors

def inverse_power_method(mat, N = N_MAX, return_eigenvector = False, shifted = True):
    '''
    Finds and returns the smallest eigenvalue of a square matrix
    '''

    np.random.seed(0)
    # x0 = np.array([row[0] for row in mat])
    x0 = np.random.random(len(mat)) + 1j * np.random.random(len(mat))
    y_n = [x0]

    for _ in range(N):
        if shifted: new_mat = mat - np.identity(len(mat)) * _rayleigh_quotient(mat, y_n[-1])
        else: new_mat = mat - np.identity(len(mat)) * _rayleigh_quotient(mat, y_n[-1])
        new_y = solve_linear_system(new_mat, y_n[-1])[0]
        new_y_normalized = new_y / np.linalg.norm(new_y)
        y_n.append(new_y_normalized)
    
    approx = _rayleigh_quotient(mat, y_n[-1])
    if return_eigenvector: return approx, y_n[-1]
    else: return approx


if __name__ == "__main__":

    mat = np.array([[5, -1, 3, 0], [-8, 26, 1903, 666], [2, 0, 7, -26], [0, 4, 4, -44]])
    print(mat)

    eigenvalues, eigenvectors = power_method_all(mat)

    print(eigenvalues)
    print(np.linalg.eigvals(mat))
