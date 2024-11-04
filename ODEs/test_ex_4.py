import numpy as np
import matplotlib.pyplot as plt
from my_odelib import unpack_V2
from good_exercise_4 import simulation
import sys

H0 = .005

def conf_0():
    # STARTING_VALS ZERO
    m0 = .3
    m1 = .3
    m2 = .3

    MASSES = np.array([m0, m1, m2])

    r0_0 = np.array([1, 0, 0], dtype = np.float128)
    r1_0 = np.array([-1, 0, 0], dtype = np.float128)
    r2_0 = np.array([0, 0, 0], dtype = np.float128)

    v0_0 = np.array([0, .15, -.15], dtype = np.float128)
    v1_0 = np.array([0, -.15, .15], dtype = np.float128)
    v2_0 = np.array([0, 0, 0], dtype = np.float128)

    START_VALS = [0., np.array([r0_0, r1_0, r2_0, v0_0, v1_0, v2_0])] # Starting values for t, x, y, z, x_dot, y_dot, z_dot

    simulation(start_vals = START_VALS, masses = MASSES, steps_number = 20_000, h = H0)


def conf_1():
    # STARTING_VALS ONE
    m0 = 1.6
    m1 = .4
    m2 = .4

    MASSES = np.array([m0, m1, m2])

    r0_0 = np.array([1, 0, 0], dtype = np.float128)
    r1_0 = np.array([-1, 0, 0], dtype = np.float128)
    r2_0 = np.array([0, 0, 0], dtype = np.float128)

    v0_0 = np.array([0, .4, 0], dtype = np.float128)
    v1_0 = np.array([0, -.8, .7], dtype = np.float128)
    v2_0 = np.array([0, -.8, -.7], dtype = np.float128)

    START_VALS = [0., np.array([r0_0, r1_0, r2_0, v0_0, v1_0, v2_0])] # Starting values for t, x, y, z, x_dot, y_dot, z_dot

    simulation(start_vals = START_VALS, masses = MASSES, steps_number = 31_600, h = H0)


def conf_2():
    # STARTING_VALS TWO
    m0 = 1.6
    m1 = .4
    m2 = .6

    MASSES = np.array([m0, m1, m2])

    r0_0 = np.array([1, 0, 0], dtype = np.float128)
    r1_0 = np.array([-1, 0, 0], dtype = np.float128)
    r2_0 = np.array([0, 0, 0], dtype = np.float128)

    v0_0 = np.array([0, .4, 0], dtype = np.float128)
    v1_0 = np.array([0, -.8, .7], dtype = np.float128)
    v2_0 = np.array([0, -.8, -.7], dtype = np.float128)

    START_VALS = [0., np.array([r0_0, r1_0, r2_0, v0_0, v1_0, v2_0])] # Starting values for t, x, y, z, x_dot, y_dot, z_dot

    simulation(start_vals = START_VALS, masses = MASSES, steps_number = 40_000, h = H0)



if __name__ == "__main__":
    '''
    Execute ex. 4 with one of many configurations (0, ...N)
    '''
    
    configs = [conf_0, conf_1, conf_2]

    if len(sys.argv) > 1:
        try: index = int(sys.argv[1])
        except (TypeError, ValueError):
            print("Bad argument")
        else: configs[index]()
    else: raise ValueError("No argument")

