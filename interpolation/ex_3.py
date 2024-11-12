import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from matrices.matrix_utils import solve_linear_system
from ex_2 import chebyshev_nodes
from first_order_spline import linear_spline
from second_order_spline import quadratic_spline

from Polynomial_classes import Polynomial
from my_interplib import *

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

INTERVAL = (-1, 1)
DELTA = 11
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def interp_and_plot(delta, spline, index, axes, func = runge, interval = INTERVAL):

    sample = np.linspace(*INTERVAL, num = delta)
    f = [func(x) for x in sample]
    arrx, arry, interp = spline(sample, f, delta, interval = interval)

    axes[index].scatter(sample, f, color = (.1, .1, .1), marker = "x", label = "Samples")
    axes[index].plot(arrx, arry, color = colors[index], alpha = 1, label = f"{spline.__name__} with {len(sample)} samples")
    axes[index].legend()

if __name__ == "__main__":
    
    fig, axes = plt.subplots(2)
    # plt.ylim(-1, 2)
    fig.set_size_inches(20, 10, forward = True)

    interp_and_plot(20, linear_spline, 0, axes)
    interp_and_plot(21, linear_spline, 1, axes)


    arrx = np.linspace(*INTERVAL, 1_000)
    for ax in axes:
        ax.plot(arrx, runge(arrx), c = (.3, .3, .3), alpha = .1, linewidth = 10, label = "Runge function")

    plt.legend()
    plt.savefig(f"interp_graphs/ex3_linear")
    plt.show()



    fig, axes = plt.subplots(2)
    # plt.ylim(-1, 2)
    fig.set_size_inches(20, 10, forward = True)

    interp_and_plot(20, quadratic_spline, 0, axes)
    interp_and_plot(21, quadratic_spline, 1, axes)


    arrx = np.linspace(*INTERVAL, 1_000)
    for ax in axes:
        ax.plot(arrx, runge(arrx), c = (.3, .3, .3), alpha = .1, linewidth = 10, label = "Runge function")

    plt.legend()
    plt.savefig(f"interp_graphs/ex3_quadratic")
    plt.show()