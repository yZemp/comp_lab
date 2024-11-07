import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("~/Documents/Programming/comp_lab/")
from interp_direct_monomial import interp_simple
from ex_1 import interp_newton
from my_interplib import *

#######################################################################
# VARS

INTERVAL = (-1, 1)
SAMPLES = 30
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def main(sample, label):
    '''
    Perform interpolation of the runge function, given a specific sample
    '''

    f = [runge(x) for x in sample]

    simple = interp_simple(sample, f)
    newton = interp_newton(sample, f)

    arrx = np.linspace(*INTERVAL, num = 5_000)

    data_simple = [simple.evaluate(el) for el in arrx]
    data_newton = [newton.evaluate(el) for el in arrx]

    fig, ax = plt.subplots(1)

    ax.scatter(sample, f, color = (.1, .1, .1), label = "Data sample", s = 70, marker = "x")
    ax.plot(arrx, runge(arrx), c = colors[0], label = "Runge function")
    ax.plot(arrx, data_simple, c = colors[1], label = "Simple interpolation")
    ax.plot(arrx, data_newton, c = colors[2], label = "Newton interpolation")

    plt.legend()
    plt.ylim(-1, 2)
    fig.suptitle(label.title() + " sampling")
    fig.set_size_inches(20, 10, forward = True)
    plt.savefig(f"interp_graphs/ex2_{label}")
    plt.show()


if __name__ == "__main__":
    # Uniform sample
    uniform = np.linspace(*INTERVAL, num = SAMPLES)
    main(uniform, "uniform")

    # Chebyshev sample
    chebyshev = chebyshev_nodes(SAMPLES)
    main(chebyshev, "chebyshev")