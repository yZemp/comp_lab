import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from my_odelib import unpack_V2

######################################################################################
# Useful vars

H0 = .005 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 20_000 # Default steps number
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def _frac(ri, rj):
    return (ri - rj) / np.power(np.linalg.norm(ri - rj), 3)

def _chi(i, R, m):
    return np.sum(np.array([- m[j] * _frac(R[i], R[j]) for j in [0, 1, 2] if j != i]), axis = 0)

def f(t, Y, m):
    x, y, z, vx, vy, vz = Y
    R = [x, y, z]
    return np.array([vx, vy, vz, _chi(0, R, m = m), _chi(1, R, m = m), _chi(2, R, m = m)])

def rk4_step(t, Y, h, known, m):
    k1 = known(t, Y, m = m)
    k2 = known(t + .5 * h, Y + .5 * h * k1, m = m)
    k3 = known(t + .5 * h, Y + .5 * h * k2, m = m)
    k4 = known(t + h, Y + h * k3, m = m)
    return (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

def solve(Y0, func, m, N, final_time = FINAL_TIME, step = rk4_step, h = H0):
    steps = [Y0]
    for _ in range(N):
        t, Y = steps[-1][0], steps[-1][1]
        new_step = [t + h, Y + h * step(t, Y, h, known = func, m = m)]
        steps.append(new_step)
    return steps


def simulation(start_vals, masses, steps_number = STEPS_NUMBER, h = H0, speed = 100):

    fig = plt.figure(figsize = (20, 20))
    ax = fig.add_subplot(111, projection = "3d")

    ax.set_xlim((-2, 2))
    ax.set_ylim((-2, 2))
    ax.set_zlim((-2, 2))

    line0, = ax.plot([], [], [], color = colors[0], label = "Body 1", lw = 2)
    line1, = ax.plot([], [], [], color = colors[1], label = "Body 2", lw = 2)
    line2, = ax.plot([], [], [], color = colors[2], label = "Body 3", lw = 2)
    
    coords = unpack_V2(solve(start_vals, func = f, m = masses, N = steps_number))
    r0 = np.array([coords[1][i] for i in range(len(coords[1]))])
    r1 = np.array([coords[2][i] for i in range(len(coords[2]))])
    r2 = np.array([coords[3][i] for i in range(len(coords[3]))])

    def update(frame):
        line0.set_data(r0[:frame, 0], r0[:frame, 1])
        line0.set_3d_properties(r0[:frame, 2])

        line1.set_data(r1[:frame, 0], r1[:frame, 1])
        line1.set_3d_properties(r1[:frame, 2])

        line2.set_data(r2[:frame, 0], r2[:frame, 1])
        line2.set_3d_properties(r2[:frame, 2])
        
        return line0, line1, line2

    ani = FuncAnimation(fig, update, frames = range(0, len(r0), speed), blit = True)

    plt.legend()
    plt.show()
