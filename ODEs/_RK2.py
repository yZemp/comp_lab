import numpy as np
import matplotlib.pyplot as plt

final_time = 30

def RK2(t0, theta0, theta_dot0, h = .001, final_time = final_time):
    arr_t = [t0]
    arr_theta = [theta0]
    arr_theta_dot = [theta_dot0]

    domain = (t0, final_time)

    for i in range(int((domain[1] - domain[0]) / h)):
        new_t = arr_t[-1] + h

        new_theta_dot = arr_theta_dot[-1] + h * f(arr_t[-1] + .5 * h, arr_theta[-1] + .5 * h * f(arr_t[-1], arr_theta[-1]))
        new_theta = arr_theta[-1] + h * arr_theta_dot[-1]

        arr_t.append(new_t)
        arr_theta.append(new_theta)
        arr_theta_dot.append(new_theta_dot)

    return [arr_t, arr_theta, arr_theta_dot]

# Analytic solution

def solution(t): return np.sin(t)


def f(t, theta): return - theta


if __name__ == "__main__":
    '''
    Solving the system Theta_dot_dot = - Theta
    with initial value (Theta(0), Theta_dot(0)) = (0, 1)
    
    This is decomposed in:
        Theta_dot(t) = phi
        phi_dot(t) = - theta
    '''

    # RUNGE-KUTTA 2 (MIDPOINT METHOD)

    solved = RK2(0., 0., 1.)
    # print(len(solved))

    fig, ax = plt.subplots(2)
    fig.set_size_inches(200, 100, forward = True)
    fig.suptitle("RK2", fontsize = 20)

    lnsp = np.linspace(0, final_time, 1000)
    ax[0].plot(lnsp, solution(lnsp), c = (0., .1, .3), label = "Correct $\Theta(t)$")
    plt.xlabel("t")

    arrx = []
    arry = []
    for x, y, _ in zip(*solved):
        arrx.append(x)
        arry.append(y)
    
    ax[0].plot(arrx, arry, c = (.1, .7, .3), label = "Euler h = .001")
    ax[0].legend()
    plt.show()
    fig.savefig("RK2_harmonic_oscillator.png", dpi = fig.dpi, bbox_inches = "tight")