import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2
import random

######################################################################################
# Useful vars

H0 = .001 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 50_000 # Default steps number
START_VALS = [0., np.array([65., 124.], dtype = np.float128)] # Starting values for t, x, y
A = 1.5
B = 1.3
C = 0.05
D = 0.07
prec = 4
epsilon = .6

def f(t, Y):
    '''
    Known term of the ODE system (in terms of t, x, y)
    args is expected to be (x, y)
    '''

    x, y = Y

    f1 = A * x - D * x * y
    f2 = - B * y + C * x * y

    return np.array([f1, f2])

def f_pesca(t, Y):
    '''
    Known term of the ODE system (in terms of t, x, y)
    args is expected to be (x, y)
    '''

    x, y = Y

    f1 = A * x - D * x * y - epsilon * x
    f2 = - B * y + C * x * y - epsilon * y

    return np.array([f1, f2])


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
    Solve a Lotka-Volterra system
    '''

    arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)

    # funcs_step = [euler_step, rk2_step, rk4_step]
    funcs_step = [rk4_step]
    colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

    for i, fun_step in enumerate(funcs_step):

        fig, ax = plt.subplots(1)
        fig.set_size_inches(20, 10, forward = True)

        coords = unpack_V2(solve(START_VALS, func = f, step = fun_step))
        # print(coords)

        ax.plot(coords[0], coords[1], c = colors[0], label = "Prede")
        ax.plot(coords[0], coords[2], c = colors[1], label = "Predatori")

        plt.legend()
        plt.xlabel("Time")
        plt.suptitle(f"{fun_step.__name__}")
        plt.savefig(f"ex_3_graphs/LV_{fun_step.__name__}.png")
        plt.show()

    fig, ax = plt.subplots(1)
    
    for x in [[x, y] for x in [2, 5] for y in np.linspace(10, 200, num = prec)]:
        print(x)
        coords = unpack_V2(solve([0., x], func = f, step = fun_step))
        ax.plot(coords[1], coords[2], c = colors[0])

    # plt.xlabel("Prede")
    # plt.ylabel("Predatori")
    # plt.suptitle("Phase plane")
    # plt.savefig(f"ex_3_graphs/LV_phaseplane_{fun_step.__name__}.png")
    # plt.show()

    
    # PESCA
    # fig, ax = plt.subplots(1)
    
    for x in [[x, y] for x in [2, 5] for y in np.linspace(10, 200, num = prec)]:
        print(x)
        coords = unpack_V2(solve([0., x], func = f_pesca, step = fun_step))
        ax.plot(coords[1], coords[2], c = colors[1])

    plt.xlabel("Prede")
    plt.ylabel("Predatori")
    plt.suptitle("Phase plane")
    plt.show()
