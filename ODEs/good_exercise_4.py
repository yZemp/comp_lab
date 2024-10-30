import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2

######################################################################################
# Useful vars

H0 = .005 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 10_000 # Default steps number
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]

# # STARTING_VALS ONE
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


# STARTING_VALS TWO
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

def _cycle(i, n = 2):
    return i if i < n else i - (i // (n + 1)) * (n + 1) 

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

    coords = unpack_V2(solve(START_VALS, func = f))
    '''
    coords should be:
        [arr_t, arr_r0, arr_r1, arr_r2, arr_v0, arr_v1, arr_v2]
        where
            arr_t is array of scalars
            others are array of 3-dim array (vectors coordinates)
    '''  

    #######################################################################################
    # First visualization
    # for j in [0, 1, 2]:
    #     fig, ax = plt.subplots(3)
    #     fig.set_size_inches(20, 10, forward = True)
        
    #     for i in [0, 1, 2]:
    #         ax[i].set_title(f"{i}-nth coord as a function of t")
    #         ax[i].plot([el for el in coords[0]], [el[i] for el in coords[j + 1]], c = colors[j], label = f"Body {j}")
    #         ax[i].legend()
    #     plt.show()
        

    #######################################################################################
    # Second visualization
    r0x = [el[0] for el in coords[1]]
    r0y = [el[1] for el in coords[1]]
    r0z = [el[2] for el in coords[1]]

    r1x = [el[0] for el in coords[2]]
    r1y = [el[1] for el in coords[2]]
    r1z = [el[2] for el in coords[2]]

    r2x = [el[0] for el in coords[3]]
    r2y = [el[1] for el in coords[3]]
    r2z = [el[2] for el in coords[3]]
    
    r0 = [r0x, r0y, r0z]
    r1 = [r1x, r1y, r1z]
    r2 = [r2x, r2y, r2z]

    r = [r0, r1, r2]

    for j in [0, 1, 2]:
        fig, ax = plt.subplots()
        fig.set_size_inches(20, 10, forward = True)
        
        plt.title(f"Evolution of motion on plane {j + 1}")
        for i, r_ in enumerate(r):
            plt.plot(r_[_cycle(j)], r_[_cycle(j + 1)], c = colors[i], label = f"Body {i}")
            
            plt.legend()
        plt.show()

    fig = plt.figure()
    ax = plt.axes(projection='3d')


    #######################################################################################
    # Third visualization    
    for i, r_ in enumerate(r):
        x, y, z = r_[0], r_[1], r_[2]
        ax.plot3D(x, y, z, c = colors[i])
    plt.show()
            

    # plt.suptitle(f"{fun_step.__name__}")
    # plt.savefig(f"ex_3_graphs/ex3_{fun_step.__name__}.png")
