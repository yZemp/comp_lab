import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 50 # Default length of approximation
STEPS_NUMBER = 5_000 # Default steps number
START_VALS = [0., np.array([0., 1.])] # Starting values for t, theta, phi
gamma = .3
A = 1.2


# Analytical solution
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


def solve(Y0, func, N = STEPS_NUMBER, final_time = FINAL_TIME, step = euler_step):
    
    steps = [Y0]

    h = final_time / N

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

    NOTE: f1, f2, f3 does NOT represent different known part of a differential system, but three different exercices
    '''


    arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)
    # plt.plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution") # I DONT HAVE IT :(

    funcs = [f1, f2, f3]


    fig, ax = plt.subplots(len(funcs))
    fig.set_size_inches(20, 10, forward = True)

    for i, fun in enumerate(funcs):
        euler_coords = unpack_V2(solve(START_VALS, func = fun))
        ax[i].plot(euler_coords[0], euler_coords[1], label = f"Euler solution for {f1.__name__}")
        ax[i].legend()
    
    plt.show()