import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse

from schr_1D import schr_solve
from Potentials import Potential_simple
from interpolation.nth_order_spline import nth_spline

#######################################################################
# VARS

N_MAX = 10_000
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]
COLOR_SHIFT_SCALE = 1
SEPARATOR = .2


A_DEF = 0.
B_DEF = 0.
N_DEF = 64
L_DEF = 1
M_DEF = 8 / L_DEF
DX_DEF = L_DEF / N_DEF
V0_DEF = -10.


def _colormap(i, norm, third = .3):
    return (min(i * COLOR_SHIFT_SCALE / norm, 1), min(i * COLOR_SHIFT_SCALE / norm, 1), third)



def ex_1(potential, label = "", number = 4):

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
    plt.savefig(f"schr_graphs/eigenstates_pdf_{label}.png")
    plt.show()


def ex_2(potential, label = "", number = 4):

    eigenvals, eigenvecs = schr_solve(N, potential, A, B, m, L, dx, N_MAX)

    arri = [i / N for i in range(N)]

    plt.figure(figsize = (20, 10))

    for i, arr in enumerate(eigenvecs):
        if i > number: continue
        plt.scatter(arri, arr, color = _colormap(i, number), alpha = .8)
        x, y, _ = nth_spline(arri, arr, N, order = 3, interval = (0, arri[-1]))
        plt.plot(x, y, color = _colormap(i, number), label = f"Eigenstate {eigenvals[i]}")

    plt.title("First eigenstates of the system")

    plt.legend()
    plt.savefig(f"schr_graphs/eigenstates_{label}.png")
    plt.show()


def ex_4(potentials, label = ""):

    arri = [i / N for i in range(N)]
    plt.figure(figsize = (20, 10))

    for i, potential in enumerate(potentials):
        eigenvals, eigenvecs = schr_solve(N, potential, A, B, m, L, dx, 2_000)

        plt.scatter(arri, eigenvecs[0], color = colors[i], alpha = .8)
        x, y, _ = nth_spline(arri, eigenvecs[0], N, order = 3, interval = (0, arri[-1]))
        plt.plot(x, y, color = colors[i], label = f"Eigenstate {eigenvals[0]}")

    plt.title("First eigenstates of the system (varying potential)")

    plt.legend()
    plt.savefig(f"schr_graphs/eigenstates_{label}.png")
    plt.show()


def ex_5(potential, label = ""):

    fig, axes = plt.subplots(3)

    for i, N in enumerate([16, 32, 64]):
        arri = [i / N for i in range(N)]
        eigenvals, eigenvecs = schr_solve(N, potential, A, B, m, L, dx, 2_000)

        for j in [0, 3, 7]:
            axes[i].scatter(arri, eigenvecs[j], color = colors[i], alpha = .8)
            x, y, _ = nth_spline(arri, eigenvecs[j], N, order = 3, interval = (0, arri[-1]))
            axes[i].plot(x, y, color = colors[i], label = f"Eigenstate {eigenvals[j]}")

        axes[i].legend()
    
    fig.set_size_inches((20, 10), forward = True)

    plt.suptitle("Eigenstates varying lattice spacing")

    plt.savefig(f"schr_graphs/eigenstates_{label}.png")
    plt.show()


if __name__ == "__main__":
    global N, L, m, dx, A, B
    
    if len(sys.argv) < 5: print("Using some default values for constants")

    '''
    Example usage: 
        python3 schr_ex.py N 64 V0 -10 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', type=int, default=N_DEF, help="Lattices number")
    parser.add_argument('--V0', type=float, default=V0_DEF, help="Potential well's depth")
    parser.add_argument('--A', type=float, default=A_DEF,)
    parser.add_argument('--B', type=float, default=B_DEF,)
    parser.add_argument('--L', type=float, default=L_DEF,)
    
    args = parser.parse_args()
    
    N = args.N
    A = args.A
    B = args.B
    L = args.L
    m = 8 / L
    dx = L / N

    print("Initiated parameters:\n", N, A, B, L, m, dx)


    potential = Potential_simple([L / 4, (3 / 4) * L], V0 = args.V0)
    inv_potential = Potential_simple([L / 4, (3 / 4) * L], V0 = - args.V0)
    potentials = [Potential_simple([L / 4, (3 / 4) * L], V0 = el) for el in [10, 40, 80]]

    # ex_1(potential, "well")
    # ex_2(potential, "well")
    # ex_2(inv_potential, "wall")
    # ex_4(potentials, "various")
    ex_5(potential, "lattice")
