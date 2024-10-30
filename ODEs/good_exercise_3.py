import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2
import random

######################################################################################
# Useful vars

H0 = .002 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 10_000 # Default steps number
START_VALS = [0., np.array([7., 12., 15.])] # Starting values for t, x, y, z

def f(t, Y):
    '''
    Known term of the ODE system (in terms of t, x, y, z)
    args is expected to be (x, y, z)
    '''

    x, y, z = Y

    f1 = 10 * (y - x)
    f2 = - x * z + 28 * x - y
    f3 = x * y - (8 / 3) * z

    return np.array([f1, f2, f3])


def euler_step(t, Y, h, known):
    return known(t, Y)

def rk2_step(t, Y, h, known):
    k1 = known(t, Y)
    k2 = known(t + .5 * h, Y + .5 * h * k1)
    return k2

def rk4_step(t, Y, h, known):
    k1 = known(t, Y)
    k2 = known(t + .5 * h, Y + .5 * h * k1)
    k3 = known(t + .5 * h, Y + .5 * h * k2)
    k4 = known(t + h, Y + h * k3)
    return (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def solve(Y0, func, N = STEPS_NUMBER, final_time = FINAL_TIME, step = euler_step, h = H0):
    
    steps = [Y0]

    # h = final_time / N

    for _ in range(N):
        t, Y = steps[-1][0], steps[-1][1]
        new_step = [t + h, Y + h * step(t, Y, h, known = func)]
        steps.append(new_step)

    return steps


if __name__ == "__main__":
    '''
    Solve a system of first order ODE
    '''


    arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)

    funcs_step = [euler_step, rk2_step, rk4_step]
    colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

    for i, fun_step in enumerate(funcs_step):

        fig, ax = plt.subplots(3)
        fig.set_size_inches(20, 10, forward = True)

        coords = unpack_V2(solve(START_VALS, func = f, step = euler_step))
        # print(coords)

        ax[0].plot(coords[1], coords[2], c = colors[i])
        ax[1].plot(coords[1], coords[3], c = colors[i])
        ax[2].plot(coords[2], coords[3], c = colors[i])

        plt.suptitle(f"{fun_step.__name__}")
        # plt.savefig(f"ex_3_graphs/ex3_{fun_step.__name__}.png")
        plt.show()

        if i == 2:
            fig = plt.figure()
            ax = plt.axes(projection = "3d")
            
            x, y, z = coords[1], coords[2], coords[3]
            ax.plot3D(x, y, z, c = colors[i])
            plt.show()
