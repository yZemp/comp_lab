import numpy as np
import matplotlib.pyplot as plt
import random as rand
import scipy.stats as sts

import sys
sys.path.append("~/Documents/Programming/comp_lab/")
# from matrices.matrix_utils import get_inverse
from Polynomial_classes import Newton_polynomial, Newton_interpolator
from interp_direct_monomial import interp_simple

# Direct approach using Newton polynomials

def _divided_diff(ei, ef, xi, xf):
    return (ef - ei) / (xf - xi)

def _arr_divided_diff(arr):
    '''
    Should receive array of array of f_i
    Returns array of divided differences
    '''

    turbomat = [arr]

    for k in range(len(arr) - 1):

        # Here x is global data
        new = np.array([_divided_diff(turbomat[k][i], turbomat[k][i + 1], x[i], x[i + k + 1]) for i in range(0, len(turbomat[k]) - 1)], np.float64)
        turbomat.append(new)

    return [row[0] for row in turbomat]


# Init data
glob_x = np.array([0, 10, 15, 20, 22.5, 30]) # i = 0, 1, 2, 3, 4, 5
glob_f = np.array([0, 227.04, 362.78, 517.35, 602.97, 901.67])

# Splitting data:
start = 0
end = 6
x = glob_x[start:end]
f = glob_f[start:end]


if __name__ == "__main__":

    f_j = _arr_divided_diff(f)

    # Building newton polynomial array ([n_0, n_1, n_2, n_3, ..., n_i])
    newton_poly_arr = [Newton_polynomial(x[0:i]) for i in range(0, len(x))]

    # Building interpolator
    print([i for i in range(start, end)], "\n", f_j, newton_poly_arr)
    newton_interpolator = Newton_interpolator(f_j, newton_poly_arr)

    arrx = np.linspace(-5, 35)

    fig, ax = plt.subplots(1)

    ax.scatter(glob_x, glob_f, color = (.9, .1, .1), s = 70, label = "Fixed points", marker = "x")

    interpolated = [newton_interpolator.evaluate(el) for el in arrx]
    ax.plot(arrx, interpolated, color = (.5, .3, .9), label = "Newton fit")


    # CONFRONTING WITH SIMPLEST FIT
    lnsp, y = interp_simple(x, f)
    # ax.plot(lnsp, y, c = (.7, .7, .7), label = "Simple fit")

    plt.legend()
    # plt.ylim(-500, 2_000)
    fig.suptitle(f"Imposing points: [{start}, {end - 1}]")
    fig.set_size_inches(20, 10, forward = True)
    plt.savefig(f"interp_graphs/interp_{start}_{end - 1}")
    plt.show()
    