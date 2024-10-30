import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2
import random

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 300 # Default length of approximation
STEPS_NUMBER = 5_000 # Default steps number
START_VALS = [0., np.array([0., 1.])] # Starting values for t, theta, phi
gamma = .2
A = 1.8

def f1(t, Y):
    '''
    Known term of the ODE (in terms of t, theta, phi)
    args is expected to be (theta, phi)
    '''

    return np.array([Y[1], - np.sin(Y[0])])


def f2(t, Y):
    '''
    Known term of the ODE (in terms of t, theta, phi)
    args is expected to be (theta, phi)
    '''

    return np.array([Y[1], - np.sin(Y[0]) - gamma * Y[1]])


def f3(t, Y):
    '''
    Known term of the ODE (in terms of t, theta, phi)
    args is expected to be (theta, phi)
    '''

    return np.array([Y[1], - np.sin(Y[0]) - gamma * Y[1] + A * np.sin((2 / 3) * t)])


def rk4_step(t, Y, h, known):
    k1 = known(t, Y)
    k2 = known(t + .5 * h, Y + .5 * h * k1)
    k3 = known(t + .5 * h, Y + .5 * h * k2)
    k4 = known(t + h, Y + h * k3)
    return (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def solve(Y0, func, N = STEPS_NUMBER, final_time = FINAL_TIME, step = rk4_step):
    
    steps = [Y0]

    # h = final_time / N
    h = H0

    for _ in range(N):
        t, Y = steps[-1][0], steps[-1][1]
        new_step = [t + h, Y + h * step(t, Y, h, known = func)]
        steps.append(new_step)

    return steps

if __name__ == "__main__":
    '''
    Solve a second order system:
        theta_dot_dot = f(t, theta)

    That can be rewritten in
        theta_dot = phi
        phi_dot = f(t, theta)

    For every step, Y is the vector (t, theta, phi) with
        t = time
        theta = height of point
        phi = first derivative of theta

    NOTE: f1, f2, f3 does NOT represent different known parts of a differential system, but three different exercices
    '''


    arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)

    funcs = [f1, f2, f3]


    for i, fun in enumerate(funcs):
        fig, ax = plt.subplots(3)
        fig.set_size_inches(20, 10, forward = True)
        
        coords = unpack_V2(solve(START_VALS, func = fun, step = rk4_step))
        
        ax[0].plot(coords[0], coords[1], c = (.5, .1, .8))
        ax[1].plot(coords[0], coords[2], c = (.5, .1, .8))
        ax[2].plot(coords[2], coords[1], c = (.5, .1, .8))

        plt.suptitle(f"{fun.__name__}")
        plt.savefig(f"ex_2_graphs/ex2_{fun.__name__}.png")
        plt.show()
