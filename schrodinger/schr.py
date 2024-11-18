import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from eigensystems.eigensystem import power_method, power_method_all, inverse_power_method
from interpolation.second_order_spline import quadratic_spline

np.set_printoptions(linewidth=np.inf)


#######################################################################
# VARS

N_MAX = 300
COLOR_SHIFT_SCALE = 10
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

A = 0.
B = 0.
N = 64
L = 1
m = 8 / L
dx = L / N
V0 = 10.

def _cround(z, threshold = 10e-8):
    if abs(z.imag) < threshold:
        return z.real
    return z

def POTENTIAL(x):
    return - V0 if L / 4. < x < (3. / 4.) * L else 0.

def init_K(n):
    K = np.diag([2] * n) + np.diag([- 1] * (n - 1), 1) + np.diag([- 1] * (n - 1), - 1)
    K[0][-1] = A
    K[-1][0] = B
    return K 

def init_V(n, potential = POTENTIAL):
    V = np.zeros((n, n))
    for i in range(0, len(V)):
        V[i][i] = potential(i * dx)
    return V


if __name__ == "__main__":
    K = init_K(N)
    V = init_V(N)
    print(K, "\n\n", V)

    bigmat = - (np.power(N, 2) / (2 * m * L)) * K + V

    eigenvals, eigenvecs = power_method_all(bigmat, N_MAX)
    eigenvals = [_cround(z) for z in eigenvals]

    print(len(eigenvals), "\n", eigenvals)

    arri = [i / N for i in range(N)]

    plt.figure(figsize = (20, 10))

    for i, arr in enumerate(eigenvecs):
        print(arri, np.abs(arr))
        if i > 5: break
        print("Autofunzioni:\n", arr)
        # plt.scatter(arri, arr, color = (i * COLOR_SHIFT_SCALE / len(eigenvecs), i * COLOR_SHIFT_SCALE / len(eigenvecs), .5), label = "")
        plt.plot(arri, arr, color = (i * COLOR_SHIFT_SCALE / len(eigenvecs), i * COLOR_SHIFT_SCALE / len(eigenvecs), .5), marker = "o", label = f"{i}")
        # x, y, _ = quadratic_spline(arri, arr, N - 1)
        # plt.plot(x, y, color = (i / 5, i / 5, .2), label = "")

    arrx = np.linspace(0, L, num = 10_000)
    arrv = [POTENTIAL(x) for x in arrx]
    # plt.plot(arrx, arrv, color = (.1, .1, .1), alpha = 1, linestyle = "dotted", label = "Potential")

    plt.title("First eigenstates of the system (Real part only)")

    plt.legend()
    plt.savefig("schr_graphs/eigenstates.png")
    plt.show()