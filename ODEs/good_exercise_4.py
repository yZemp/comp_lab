import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2
import random

######################################################################################
# Useful vars

H0 = .005 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 5_000 # Default steps number


# STARTING_VALS ONE
# m0 = .3
# m1 = .3
# m2 = .3

# m = np.array([m0, m1, m2])

# r0_0 = np.array([1, 0, 0], dtype = np.float128)
# r1_0 = np.array([-1, 0, 0], dtype = np.float128)
# r2_0 = np.array([0, 0, 0], dtype = np.float128)

# v0_0 = np.array([0, .15, -.15], dtype = np.float128)
# v1_0 = np.array([0, -.15, .15], dtype = np.float128)
# v2_0 = np.array([0, 0, 0], dtype = np.float128)

# START_VALS = [0., np.array([r0_0, r1_0, r2_0, v0_0, v1_0, v2_0])] # Starting values for t, x, y, z, x_dot, y_dot, z_dot


# STARTING_VALS ONE
m0 = 1.6
m1 = .4
m2 = .4

m = np.array([m0, m1, m2])

r0_0 = np.array([1, 0, 0], dtype = np.float128)
r1_0 = np.array([-1, 0, 0], dtype = np.float128)
r2_0 = np.array([0, 0, 0], dtype = np.float128)

v0_0 = np.array([0, .4, 0], dtype = np.float128)
v1_0 = np.array([0, -.8, .7], dtype = np.float128)
v2_0 = np.array([0, -.8, -.7], dtype = np.float128)

START_VALS = [0., np.array([r0_0, r1_0, r2_0, v0_0, v1_0, v2_0])] # Starting values for t, x, y, z, x_dot, y_dot, z_dot


############################################################################################3

def _frac(ri, rj):
    return (ri - rj) / np.power(np.linalg.norm(ri - rj), 3)

def _chi(i, R):
    '''
    Should be: R = [r0, r1, r2]
    '''
    return np.sum(np.array([- m[j] * _frac(R[i], R[j]) for j in [0, 1, 2] if j != i]), axis = 0)

def f(t, Y):
    '''
    Known term of the ODE system (in terms of t, x, y, z, vx, vy, vz)
    args is expected to be (x, y, z, vx, vy, vz)
    '''

    x, y, z, vx, vy, vz = Y
    R = [x, y, z]

    return np.array([vx, vy, vz, _chi(0, R), _chi(1, R), _chi(2, R)])


def rk4_step(t, Y, h, known):
    k1 = known(t, Y)
    k2 = known(t + .5 * h, Y + .5 * h * k1)
    k3 = known(t + .5 * h, Y + .5 * h * k2)
    k4 = known(t + h, Y + h * k3)
    return (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def solve(Y0, func, N = STEPS_NUMBER, final_time = FINAL_TIME, step = rk4_step, h = H0):
    
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

    colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

    fig, ax = plt.subplots()
    fig.set_size_inches(20, 10, forward = True)

    coords = unpack_V2(solve(START_VALS, func = f))
    # print(len(coords[4]))

    ax.plot(coords[1], coords[4], c = colors[0], label = "Body 1")
    ax.plot(coords[2], coords[5], c = colors[1], label = "Body 2")
    ax.plot(coords[3], coords[6], c = colors[2], label = "Body 3")

    # plt.suptitle(f"{fun_step.__name__}")
    # plt.savefig(f"ex_3_graphs/ex3_{fun_step.__name__}.png")
    ax.legend()
    plt.show()
