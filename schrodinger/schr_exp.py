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


def main(N_max = N_MAX):
    K = init_K(N)
    V = init_V(N)
    print(K, "\n\n", V)

    bigmat = - (np.power(N, 2) / (2 * m * L)) * K + V

    eigenvals, eigenvecs = power_method_all(bigmat, N_max)
    eigenvals = [_cround(z) for z in eigenvals]
    eigenvals, eigenvecs = zip(*sorted(zip(eigenvals, eigenvecs), reverse = True))

    print(len(eigenvals), "\n", eigenvals)
    # print(len(eigenvals))

    arri = [i / N for i in range(N)]

    plt.figure(figsize = (20, 10))

    for i, arr in enumerate(eigenvecs):
        # print(arri, np.abs(arr))
        # if i % 5 != 0 or i >= 10: continue
        if i > 0: continue
        # print("Autofunzioni:\n", arr)
        plt.scatter(arri, arr, color = (i * COLOR_SHIFT_SCALE / len(eigenvecs), i * COLOR_SHIFT_SCALE / len(eigenvecs), .5), alpha = .8)
        # plt.plot(arri, arr, color = (i * COLOR_SHIFT_SCALE / len(eigenvecs), i * COLOR_SHIFT_SCALE / len(eigenvecs), .5), marker = "o", label = f"{i}")
        x, y, _ = nth_spline(arri, arr, N, order = 3, interval = (0, (N - 1) * dx))
        plt.plot(x, y, color = (i * COLOR_SHIFT_SCALE / len(eigenvecs), i * COLOR_SHIFT_SCALE / len(eigenvecs), .5), label = f"Eigenstate {eigenvals[i]}")

    arrx = np.linspace(0, L, num = 10_000)
    arrv = [POTENTIAL(x) for x in arrx]
    # plt.plot(arrx, arrv, color = (.1, .1, .1), alpha = 1, linestyle = "dotted", label = "Potential")

    plt.title("First eigenstates of the system (Real part only)")

    plt.legend()
    # plt.savefig("schr_graphs/eigenstates.png")
    plt.show()


def show_convergence(N_max, ax):
    K = init_K(N)
    V = init_V(N)
    print(K, "\n\n", V)

    bigmat = - (np.power(N, 2) / (2 * m * L)) * K + V

    eigenvals, eigenvecs = power_method_all(bigmat, N_max)
    eigenvals = [_cround(z) for z in eigenvals]

    print(len(eigenvals), "\n", eigenvals)
    # print(len(eigenvals))

    arri = [i / N for i in range(N)]

    for i, arr in enumerate(eigenvecs):
        # print(arri, np.abs(arr))
        # if i % 5 != 0 or i >= 10: continue
        if i > 0: continue
        # print("Autofunzioni:\n", arr)
        ax.scatter(arri, arr, color = (N_max / 500, N_MAX / 500, .5), alpha = .8)
        # plt.plot(arri, arr, color = (i * COLOR_SHIFT_SCALE / len(eigenvecs), i * COLOR_SHIFT_SCALE / len(eigenvecs), .5), marker = "o", label = f"{i}")
        x, y, _ = nth_spline(arri, arr, N, order = 3, interval = (0, (N - 1) * dx))
        ax.plot(x, y, color = (N_max / 500, N_MAX / 500, .5), label = f"Eigenstate - {N_max} iter")



if __name__ == "__main__":
    main()


    # fig, ax = plt.subplots(1)
    # fig.set_size_inches((20, 10))

    # for n in np.linspace(1, 500, num = 10, dtype = int):
    #     show_convergence(n, ax)
    
    # plt.title("First eigenstates of the system (Real part only)")
    # plt.legend()
    # plt.show()
    