import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2
import random

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 5_000 # Default steps number
START_VALS = [0., np.array([0., 1.])] # Starting values for t, theta, phi

# Analytical solution
def solution(t): return np.sin(t)

def f(t, Y):
    '''
    Known term of the ODE (in terms of t, theta, phi)
    args is expected to be (theta, phi)
    '''

    return np.array([Y[1], -Y[0]])


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
    Solve the simple system:
        y_dot_dot = f(x, y)

    That can be rewritten in
        y_dot = z
        z_dot = f(x, y)

    For every step, Y is the vector (x, y, z) with
        x = time
        y = height of point
        z = first derivative of y
    '''

    # FIXME: No way this stuff is working correctly

    euler_coords = unpack_V2(solve(START_VALS, func = f, step = euler_step))
    rk2_coords = unpack_V2(solve(START_VALS, func = f, step = rk2_step))
    rk4_coords = unpack_V2(solve(START_VALS, func = f, step = rk4_step))


    ###############################################################################################
    # EX. 1)

    arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)

    fig, ax = plt.subplots(1)
    fig.set_size_inches(20, 10, forward = True)

    ax.plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution")


    # Euler method
    ax.plot(euler_coords[0], euler_coords[1], marker = "", c = (.1, .7, .1), label = "Euler method")

    # RK2 method
    ax.plot(rk2_coords[0], rk2_coords[1], marker = "", c = (.1, .1, .9), label = "RK2 method")

    # RK4 method
    ax.plot(rk4_coords[0], rk4_coords[1], marker = "", c = (.8, .1, .3), label = "RK4 method")


    plt.title("EX. 1.: solutions")
    plt.ylim(-2, 2)
    fig.legend()
    plt.savefig("ex_1_graphs/ex_1_1.png")
    plt.show()


    ###############################################################################################
    # EX. 2)

    arrh = np.geomspace(1, .01, endpoint = True, num = 50)

    euler_errs = []
    rk2_errs = []
    rk4_errs = []
    
    for h in arrh:    
        print(h)

        euler_coords = unpack_V2(solve(START_VALS, func = f, step = euler_step, h = h))
        rk2_coords = unpack_V2(solve(START_VALS, func = f, step = rk2_step, h = h))
        rk4_coords = unpack_V2(solve(START_VALS, func = f, step = rk4_step, h = h))

        # euler_errs.append(abs(np.array(euler_coords[1][-1]) - np.array(solution(euler_coords[0][-1]))))
        # rk2_errs.append(abs(np.array(rk2_coords[1][-1]) - np.array(solution(euler_coords[0][-1]))))
        # rk4_errs.append(abs(np.array(rk4_coords[1][-1]) - np.array(solution(euler_coords[0][-1]))))

        euler_errs.append(np.max([abs(eta - y) for eta, y in zip(euler_coords[0], euler_coords[1])]))
        rk2_errs.append(np.max([abs(eta - y) for eta, y in zip(euler_coords[0], rk2_coords[1])]))
        rk4_errs.append(np.max([abs(eta - y) for eta, y in zip(euler_coords[0], rk4_coords[1])]))

        # euler_errs.append(np.sum([abs(eta - y) for eta, y in zip(euler_coords[0], euler_coords[1])]) / h)
        # rk2_errs.append(np.sum([abs(eta - y) for eta, y in zip(euler_coords[0], rk2_coords[1])]) / h)
        # rk4_errs.append(np.sum([abs(eta - y) for eta, y in zip(euler_coords[0], rk4_coords[1])]) / h)

    fig, ax = plt.subplots(1)
    fig.set_size_inches(20, 10, forward = True)

    ax.plot(arrh, euler_errs, c = (.8, .2, .1), label = "Euler errors")
    ax.plot(arrh, rk2_errs, c = (.1, .8, .2), label = "Rk2 errors")
    ax.plot(arrh, rk4_errs, c = (.2, .1, .8), label = "Rk4 errors")

    plt.ylim(10e0, 10e23)
    plt.xlabel("h")
    # ax.gca().invert_xaxis()
    plt.title("EX. 2.: errors")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig("ex_1_graphs/ex_1_2.png")
    plt.show()
