import numpy as np
import matplotlib.pyplot as plt

final_time = 10

def euler(t0, theta0, theta_dot0, h = .01, final_time = final_time):
    arr_t = [t0]
    arr_theta = [theta0]
    arr_theta_dot = [theta_dot0]

    domain = (t0, final_time)

    for i in range(int((domain[1] - domain[0]) / h)):
        new_t = arr_t[-1] + h

        new_theta_dot = arr_theta_dot[-1] + h * f(arr_t[-1], arr_theta[-1])
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

    # FIRST EQ
    solved = euler(0., 0., 1.)
    # print(len(solved))

    lnsp = np.linspace(0, final_time, 1000)
    plt.plot(lnsp, solution(lnsp), c = (0., .1, .3), label = "Correct $\Theta(t)$")
    plt.xlabel("t")
    # plt.yscale("log")

    for x, y, _ in zip(*solved):
        plt.scatter(x, y, c = (.6, .2, .1), marker = "x")

    plt.legend()
    plt.show()