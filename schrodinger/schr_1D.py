import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys

sys.path.append("..")
from eigensystems.eigensystem import power_method_all
from Potentials import Potential_simple
from interpolation.nth_order_spline import nth_spline
from typing import Callable

np.set_printoptions(linewidth=np.inf)


#######################################################################
# VARS

N_MAX = 5_000
SEPARATOR = .2
COLOR_SHIFT_SCALE = 1

def _colormap(i, norm, third = .3):
    return (min(i * COLOR_SHIFT_SCALE / norm, 1), min(i * COLOR_SHIFT_SCALE / norm, 1), third)


def _cround(z, threshold = 10e-8):
    if abs(z.imag) < threshold:
        return z.real
    return z

def init_K(n, A, B):
    K = np.diag([- 2] * n) + np.diag([1] * (n - 1), 1) + np.diag([1] * (n - 1), - 1)
    K[0][-1] = A
    K[-1][0] = B
    return K 

def init_V(n, potential, dx):
    V = np.zeros((n, n))
    for i in range(0, n):
        V[i][i] = potential(i * dx)
    return V

def schr_solve(N: int, potential: Callable, A: float, B: float, m: float, L: float, dx: float, n_max: int = N_MAX):
    K = init_K(N, A, B)
    V = init_V(N, potential, dx)
    print(K, "\n\n", V)

    bigmat = - (np.power(N, 2) / (2 * m * L)) * K + V * L

    eigenvals, eigenvecs = power_method_all(bigmat, n_max)
    eigenvals = np.array([_cround(z) for z in eigenvals]) / L
    eigenvals, eigenvecs = map(np.array, zip(*sorted(zip(eigenvals, eigenvecs), reverse = False)))

    return eigenvals, eigenvecs


if __name__ == "__main__":
    A = 0.
    B = 0.
    L = 1.
    number = 3
    N = 64
    m = 8 / L
    dx = L / N
    potential = Potential_simple([L / 4, (3 / 4) * L], V0 = - 10)

    eigenvals, eigenvecs = schr_solve(N, potential, A, B, m, L, dx, N_MAX)
    pdfs = np.power(abs(eigenvecs), 2)

    arri = [i / N for i in range(N)]

    plt.figure(figsize = (20, 10))

    for i, arr in enumerate(pdfs):
        if i > number: continue
        plt.scatter(arri, eigenvals[i] + m * L * arr, color = _colormap(i, number), alpha = .8)
        x, y, _ = nth_spline(arri, eigenvals[i] + m * L * arr, N, order = 3, interval = (0, arri[-1]))
        plt.plot(x, y, color = _colormap(i, number), label = f"Eigenstate {eigenvals[i]}")

    arrx = np.linspace(- SEPARATOR, 1 + SEPARATOR, num = 10_000)
    plt.plot(arrx, potential(arrx), color = (.1, .1, .1), alpha = 1, linestyle = "dotted", label = "Potential")

    plt.title("First eigenstates of the system (pdf)")
    plt.legend()
    plt.show()
