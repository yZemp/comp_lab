import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys

sys.path.append("..")
from eigensystems.eigensystem import power_method_all
from typing import Callable

np.set_printoptions(linewidth=np.inf)


#######################################################################
# VARS

N_MAX = 5_000

def _cround(z, threshold = 10e-8):
    if abs(z.imag) < threshold:
        return z.real
    return z

def ndim_laplacian(x, spacings):
    pass

def init_K(n, A, B):
    K = np.diag([- 2] * n) + np.diag([1] * (n - 1), 1) + np.diag([1] * (n - 1), - 1)
    K[0][-1] = A
    K[-1][0] = B
    return K 

def init_V(n, potential, dx):
    V = np.zeros((n, n))
    for i in range(0, len(V)):
        V[i][i] = potential(i * dx)
    return V

def init(potential: Callable, A: float, B: float, m: float, L: float, spacings: list[float], n_max: int = N_MAX):
    '''
    potential: potential function to assume in schr equation
    A, B: boundary conditions
    m: mass of the particle
    L: Length of the interval (starting at 0)
    spacings: precision of approx on various axis. usually len(lat_spacing) == 1, 2, 3.
        Think about it as dx, dy, dz length
        Dimension of the problem is len(spacings) 
    n_max = max iteration of the eigensolver
    '''
    
    K = init_K(N, A, B)
    V = init_V(N, potential, dx)
    print(K, "\n\n", V)

    bigmat = - (np.power(N, 2) / (2 * m * L)) * K + V * L

    eigenvals, eigenvecs = power_method_all(bigmat, n_max)
    eigenvals = np.array([_cround(z) for z in eigenvals]) / L
    eigenvals, eigenvecs = map(np.array, zip(*sorted(zip(eigenvals, eigenvecs), reverse = False)))

    return eigenvals, eigenvecs


if __name__ == "__main__":
    pass