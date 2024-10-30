import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 20 # Default length of approximation
STEPS_NUMBER = 5_000 # Default steps number
START_VALS = [0., np.array([0., 1.])] # Starting values for t, etc

# Analytical solution
def solution(t): return np.sin(t)


# Known term of ODE (in terms of t, z0, z1, ... , zN, zN_dot)
def f(t, Y):
    '''
    This is the function that determines the last line of the system
    (See system_right)

    NOTE: this returns a scalar, and that's correct
    '''

    return - Y[0]


def system_right(t, Y, known):
    '''
    This represent the "right part" of the differential equation system of first order ODE
    known is the function that determines the last line of the system

    Known term of the ODE (in terms of t, z0, z1, ... , zN, zN_dot)
    args is expected to be z0, z1, ... , zN, zN_dot
    '''

    return np.array([*Y[1:], known(t, Y)])


def euler_step(t, Y, h, known):
    return system_right(t, Y, known)



def solve(Y0, func, N = STEPS_NUMBER, final_time = FINAL_TIME, step = euler_step):
    
    steps = [Y0]

    # h = final_time / N
    h = H0

    while True:
        t = steps[-1][0]
        Y = steps[-1][1]
        if t > final_time: return steps
        new_step = [t + h, Y + h * step(t, Y, h, known = func)]
        steps.append(new_step)



if __name__ == "__main__":
    '''

    FIXME: NOTHING WORKS HERE

    NOTE: THIS USES EULER'S METHOD ONLY

    Solve a generic ODE:
        y_(N+1)dot = f(x, y, y_dot, ... , y_(N)dot)

    That can be rewritten in
        y = z0
        y_dot = z1
        z1_dot = z2
        ...
        zN_dot = f(t, z0, z1, z2, ... , zN)

    For every step, Y is the vector (t, z0, z1, z2, ... , zN, zN_dot) with
        t = time
        z0 = height of point
        z_i = i-nth derivative of y
        zN_dot = highest order derivative of y
    '''

    arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)
    # plt.plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution")

    coords = unpack_V2(solve(START_VALS, func = f))
    plt.plot(coords[0], coords[1], c = (.1, .4, .2), label = f"Numeric solution h = {H0}")

    plt.plot(coords[0], solution(coords[0]), c = (.8, .6, .1), label = "Analytic solution")

    plt.legend()
    plt.show()