import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

from ODEs.ex_1_noclueifitworks import rk4_method, unpack

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 50 # Default length of approximation
START_VALS = (0., 0., 1.) # Starting values for x, y, z
A = 1.2
gamma = .9
color = (.1, .5, .2)

def fA(x, y):
    '''
    This represent known element of the differential equation
    '''

    return - np.sin(y)



def fC(x, y):
    '''
    This represent known element of the differential equation
    '''

    return - np.sin(y) + A * np.sin((2 * x) / 3)


def solve_and_plot(f_letter, letter, col = color):

    # Solving (using RK4)

    arrx = np.linspace(0., FINAL_TIME, 10_000)
    # plt.plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution")

    rk4_solved = rk4_method(START_VALS, h = H0, final_time = FINAL_TIME, f = f_letter)
    rk4_coords = unpack(rk4_solved)


    # Plotting 

    fig, ax = plt.subplots(3)
    fig.set_size_inches(20, 10, forward = True)
    fig.suptitle(letter)

    ax[0].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    ax[0].xaxis.set_major_locator(tck.MultipleLocator(base = 1.0))

    ax[1].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    ax[1].xaxis.set_major_locator(tck.MultipleLocator(base = 1.0))


    ax[0].plot(rk4_coords[0], rk4_coords[1], marker = "", c = col, label = "$\Theta (t)$")
    ax[0].set_xlabel("t")
    ax[0].legend()

    ax[1].plot(rk4_coords[0], rk4_coords[2], marker = "", c = col, label = "$\dot{\Theta} (t)$")
    ax[1].set_xlabel("t")
    ax[1].legend()

    ax[2].plot(rk4_coords[1], rk4_coords[2], marker = "", c = col, label = "$\dot{\Theta} (\Theta)$")
    ax[2].set_xlabel("$ \Theta $")
    ax[2].legend()


    plt.savefig(f"ex_2_graphs/{letter}.png")
    plt.close()



if __name__ == "__main__":

    solve_and_plot(fA, "A", col = (.1, .5, .2))

    solve_and_plot(fC, "C", col = (.6, .2, .1))
