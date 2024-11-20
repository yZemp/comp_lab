import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys, os

sys.path.append("..")
from eigensystems.eigensystem import power_method_all
from interpolation.nth_order_spline import nth_spline

np.set_printoptions(linewidth=np.inf)


#######################################################################
# VARS

N_MAX = 5_000
COLOR_SHIFT_SCALE = 20
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

A = 0.
B = 0.
N = 64
L = 1
m = 8 / L
dx = L / N
V0 = 10.

def _colormap(i, norm, third = .5):
    return (i * COLOR_SHIFT_SCALE / norm, i * COLOR_SHIFT_SCALE / norm, third)

def _cround(z, threshold = 10e-8):
    if abs(z.imag) < threshold:
        return z.real
    return z

def POTENTIAL(x):
    return - V0 if L / 4. < x < (3. / 4.) * L else 0.

def init_K(n):
    K = np.diag([- 2] * n) + np.diag([1] * (n - 1), 1) + np.diag([1] * (n - 1), - 1)
    K[0][-1] = A
    K[-1][0] = B
    return K 

def init_V(n, potential = POTENTIAL):
    V = np.zeros((n, n))
    for i in range(0, len(V)):
        V[i][i] = potential(i * dx)
    return V


def main(N_max = N_MAX):
    K = init_K(N)
    V = init_V(N)
    print(K, "\n\n", V)

    bigmat = - (np.power(N, 2) / (2 * m * L)) * K + V * L

    eigenvals, eigenvecs = power_method_all(bigmat, N_max)
    # eigenvals, eigenvecs = np.linalg.eig(bigmat)
    eigenvals = np.array([_cround(z) for z in eigenvals]) / L
    eigenvals, eigenvecs = map(np.array, zip(*sorted(zip(eigenvals, eigenvecs), reverse = False)))
    pdf = np.power(abs(eigenvecs), 2)

    print(len(eigenvals), "\n", eigenvals)
    # print(len(eigenvals))

    arri = [i / N for i in range(N)]

    plt.figure(figsize = (20, 10))

    for i, arr in enumerate(pdf):
        # print(arri, np.abs(arr))
        # if i % 5 != 0 or i >= 10: continue
        if i > 2: continue
        # print("Autofunzioni:\n", arr)
        plt.scatter(arri, eigenvals[i] + m * L * arr, color = _colormap(i, len(arr)), alpha = .8)
        # plt.plot(arri, arr, color = (i * COLOR_SHIFT_SCALE / len(eigenvecs), i * COLOR_SHIFT_SCALE / len(eigenvecs), .5), marker = "o", label = f"{i}")
        x, y, _ = nth_spline(arri, eigenvals[i] + m * L * arr, N, order = 3, interval = (0, (N - 1) * dx))
        plt.plot(x, y, color = _colormap(i, len(arr)), label = f"Eigenstate {eigenvals[i]}")

    arrx = np.linspace(0, L, num = 10_000)
    arrv = [POTENTIAL(x) for x in arrx]
    plt.plot(arrx, arrv, color = (.1, .1, .1), alpha = 1, linestyle = "dotted", label = "Potential")

    plt.title("First eigenstates of the system (pdf)")

    plt.legend()
    plt.savefig("schr_graphs/eigenstates.png")
    plt.show()



if __name__ == "__main__":
    main()