import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys

sys.path.append("..")
from eigensystems.eigensystem import power_method_all
from schrodinger.Potentials import Potential, Potential_simple
from typing import Callable

np.set_printoptions(linewidth=np.inf)


#######################################################################
# VARS

N_MAX = 5_000
A = 0.
B = 0.
N = 64
L = 1
m = 8 / L
dx = L / N
V0_DEFAULT = 10.
POTENTIAL = Potential_simple([0, L], V0 = - V0_DEFAULT)

def _cround(z, threshold = 10e-8):
    if abs(z.imag) < threshold:
        return z.real
    return z

def init_K(n, A, B):
    K = np.diag([- 2] * n) + np.diag([1] * (n - 1), 1) + np.diag([1] * (n - 1), - 1)
    K[0][-1] = A
    K[-1][0] = B
    return K 

def init_V(n, potential):
    V = np.zeros((n, n))
    for i in range(0, len(V)):
        V[i][i] = potential(i * dx)
    return V

def init(N: int = N, potential: Callable = POTENTIAL, A: float = A, B: float = B):
    K = init_K(N, A, B)
    V = init_V(N, potential)
    print(K, "\n\n", V)

    bigmat = - (np.power(N, 2) / (2 * m * L)) * K + V * L

    eigenvals, eigenvecs = power_method_all(bigmat, N_MAX)
    eigenvals = np.array([_cround(z) for z in eigenvals]) / L
    eigenvals, eigenvecs = map(np.array, zip(*sorted(zip(eigenvals, eigenvecs), reverse = False)))

    return eigenvals, eigenvecs


if __name__ == "__main__":
    pass