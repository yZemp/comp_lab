import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack, unpack_V2
from scipy.constants import G

######################################################################################
# Useful vars

H0 = .005 # Default precision of step
FINAL_TIME = 30 # Default length of approximation
STEPS_NUMBER = 20_000 # Default steps number
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]
E = []

############################################################################################3

def _calc_energy(Y, m):
    t = Y[0]
    x, y, z, vx, vy, vz = Y[1]
    r = [x, y, z]
    v = [vx, vy, vz]
    U = - (np.sum(np.array([(m[i] * m[j]) / np.linalg.norm(r[i] - r[j]) for i in [0, 1, 2] for j in range(i + 1, 3)])))
    K = np.sum([.5 * m_ * np.power(np.linalg.norm(v_), 2) for m_, v_ in zip(m, v)])
    return np.array([t, U, K, U + K])


def _cycle(i, n = 2):
    return i if i < n else i - (i // (n + 1)) * (n + 1) 

def _frac(ri, rj):
    return (ri - rj) / np.power(np.linalg.norm(ri - rj), 3)

def _chi(i, R, m):
    '''
    Should be: R = [r0, r1, r2]
    '''
    return np.sum(np.array([- m[j] * _frac(R[i], R[j]) for j in [0, 1, 2] if j != i]), axis = 0)

def f(t, Y, m):
    '''
    Known term of the ODE system (in terms of t, x, y, z, vx, vy, vz)
    args is expected to be (x, y, z, vx, vy, vz)
    '''

    x, y, z, vx, vy, vz = Y
    R = [x, y, z]

    return np.array([vx, vy, vz, _chi(0, R, m = m), _chi(1, R, m = m), _chi(2, R, m = m)])


def rk4_step(t, Y, h, known, m):
    k1 = known(t, Y, m = m)
    k2 = known(t + .5 * h, Y + .5 * h * k1, m = m)
    k3 = known(t + .5 * h, Y + .5 * h * k2, m = m)
    k4 = known(t + h, Y + h * k3, m = m)
    return (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def solve(Y0, m, N, func = f, final_time = FINAL_TIME, step = rk4_step, h = H0):
    
    steps = [Y0]

    # h = final_time / N

    for _ in range(N):
        t, Y = steps[-1][0], steps[-1][1]
        new_step = [t + h, Y + h * step(t, Y, h, known = func, m = m)]
        E.append(_calc_energy([new_step[0], new_step[1]], m = m))
        steps.append(new_step)

    return steps



def simulation(start_vals, masses, steps_number = STEPS_NUMBER, h = H0):
    arrx = np.linspace(0., STEPS_NUMBER * h, 10_000)

    E.append(_calc_energy(start_vals, m = masses))
    coords = unpack_V2(solve(start_vals, m = masses, N = steps_number, h = h))
    # print(E)

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
            plt.plot(r_[_cycle(j)], r_[_cycle(j + 1)], c = colors[i], label = f"Body {i + 1}")
            
            plt.legend()
        plt.show()

    #######################################################################################
    # Third visualization

    fig = plt.figure()
    ax = plt.axes(projection = "3d")

    for i, r_ in enumerate(r):
        x, y, z = r_[0], r_[1], r_[2]
        ax.plot3D(x, y, z, c = colors[i])
    plt.show()



    #######################################################################################
    # Energy visualization

    plt.figure()

    energy = unpack(E)

    plt.plot(energy[0], energy[1], c = colors[0], label = "U energy")
    plt.plot(energy[0], energy[2], c = colors[1], label = "K energy")
    plt.plot(energy[0], energy[3], c = colors[2], label = "Total energy")

    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.title("Energy Evolution over Time")
    plt.legend()
    plt.show()

    



if __name__ == "__main__":
    '''
    Solve a system of first order ODE
    '''

    '''
    coords should be:
        [arr_t, arr_r0, arr_r1, arr_r2, arr_v0, arr_v1, arr_v2]
        where
            arr_t is array of scalars
            others are array of 3-dim array (vectors coordinates)
    '''  
    pass