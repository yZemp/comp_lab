import numpy as np
import matplotlib.pyplot as plt
from nth_order import solve
from my_odelib import unpack_V2

H0 = .01 # Default precision of step
FINAL_TIME = 20 # Default length of approximation
STEPS_NUMBER = 5_000 # Default steps number
START_VALS = [0., np.array([1., 0., 0.])] # Starting values for t, etc

# Analytical solution
def solution(t): return np.exp(t) - 2 * np.exp(2 * t) + 3 * np.exp(3 * t)


# Known term of ODE (in terms of t, z1, ... , zN)
def f(t, Y):
    '''
    This is the function that determines the last line of the system
    (See system_right)

    NOTE: this returns a scalar, and that's correct
    '''

    return + 6 * Y[0] - 11 * Y[1] + 6 * Y[2]


if __name__ == "__main__":

    arrx = np.linspace(0., FINAL_TIME, 10_000)
    # plt.plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution")

    coords = unpack_V2(solve(START_VALS, func = f))
    plt.plot(coords[0], coords[1], c = (.1, .4, .2), label = f"Numeric solution h = {H0}")
    plt.ylim(-40, 200)

    plt.plot(arrx, solution(arrx), c = (.8, .6, .1), label = "Analytic solution")

    plt.legend()
    plt.show()