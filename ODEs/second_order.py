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
gamma = .3
A = 1.2


# Analytical solution
# def solution(t): return 4 * np.arctan(np.exp(t))
# def solution(t): return np.sin(t)

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


def solve(Y0, func, N = STEPS_NUMBER, final_time = FINAL_TIME, step = euler_step):
    
    steps = [Y0]

    # h = final_time / N
    h = H0

    while True:
        t, Y = steps[-1][0], steps[-1][1]
        if t > final_time: return steps
        new_step = [t + h, Y + h * step(t, Y, h, known = func)]
        steps.append(new_step)


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


    funcs_step = [euler_step, rk2_step, rk4_step]

    for funstep in funcs_step: 

        arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)

        funcs = [f1, f2, f3]


        fig, ax = plt.subplots(len(funcs))
        fig.set_size_inches(20, 10, forward = True)

        random.seed(.1)
        data = []
        for i in range(10):
            data.append([random.random() * 2, random.random() * 2])

        for i, fun in enumerate(funcs):
            for g, a in data:
                gamma = g
                A = a
                euler_coords = unpack_V2(solve(START_VALS, func = fun, step = funstep))
                ax[i].plot(euler_coords[0], euler_coords[1])
            ax[i].set_title(f"Rk4 solution for {fun.__name__}")

        # ax[0].plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution") # I DONT HAVE IT :(

        ax[0].set_ylim(-20, 20)
        ax[1].set_ylim(-.20, .20)
        ax[2].set_ylim(-5, 5)
        plt.savefig(f"ex_2_graphs/ex2_{funstep.__name__}.png")
        plt.show()