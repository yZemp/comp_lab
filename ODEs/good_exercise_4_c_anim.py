import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from my_odelib import unpack_V2

######################################################################################
# Useful vars

H0 = .005 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 40_000 # Default steps number
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

# STARTING_VALS TWO
m0 = 1.6
m1 = .4
m2 = .6

m = np.array([m0, m1, m2])

r0_0 = np.array([1, 0, 0], dtype = np.float128)
r1_0 = np.array([-1, 0, 0], dtype = np.float128)
r2_0 = np.array([0, 0, 0], dtype = np.float128)

v0_0 = np.array([0, .4, 0], dtype = np.float128)
v1_0 = np.array([0, -.8, .7], dtype = np.float128)
v2_0 = np.array([0, -.8, -.7], dtype = np.float128)

START_VALS = [0., np.array([r0_0, r1_0, r2_0, v0_0, v1_0, v2_0])] # Starting values for t, x, y, z, x_dot, y_dot, z_dot



def _frac(ri, rj):
    return (ri - rj) / np.power(np.linalg.norm(ri - rj), 3)

def _chi(i, R):
    return np.sum(np.array([- m[j] * _frac(R[i], R[j]) for j in [0, 1, 2] if j != i]), axis = 0)

def f(t, Y):
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
    for _ in range(N):
        t, Y = steps[-1][0], steps[-1][1]
        new_step = [t + h, Y + h * step(t, Y, h, known = func)]
        steps.append(new_step)
    return steps

coords = unpack_V2(solve(START_VALS, func=f))
r0 = np.array([coords[1][i] for i in range(len(coords[1]))])
r1 = np.array([coords[2][i] for i in range(len(coords[2]))])
r2 = np.array([coords[3][i] for i in range(len(coords[3]))])

fig = plt.figure(figsize = (20, 20))
ax = fig.add_subplot(111, projection = "3d")

ax.set_xlim((-2, 2))
ax.set_ylim((-2, 2))
ax.set_zlim((-2, 2))

line0, = ax.plot([], [], [], color = colors[0], label = "Body 1", lw = 2)
line1, = ax.plot([], [], [], color = colors[1], label = "Body 2", lw = 2)
line2, = ax.plot([], [], [], color = colors[2], label = "Body 3", lw = 2)

def update(frame):
    line0.set_data(r0[:frame, 0], r0[:frame, 1])
    line0.set_3d_properties(r0[:frame, 2])

    line1.set_data(r1[:frame, 0], r1[:frame, 1])
    line1.set_3d_properties(r1[:frame, 2])

    line2.set_data(r2[:frame, 0], r2[:frame, 1])
    line2.set_3d_properties(r2[:frame, 2])
    
    return line0, line1, line2

ani = FuncAnimation(fig, update, frames = range(0, len(r0), 100), blit = True)

plt.legend()
plt.show()
