import numpy as np
import matplotlib.pyplot as plt

final_time = 10

def euler(init_value: tuple, h = .1, final_time = final_time):
    solved = [init_value]
    domain = (init_value[0], final_time)

    for i in range(int((domain[1] - domain[0]) / h)):
        new_x = solved[-1][0] + h
        new_y = solved[-1][1] + h * f(solved[-1][0], solved[-1][1])
        solved.append((new_x, new_y))

    return solved


# Analytic solution

def solution(t): return np.exp(k * t)

k = .7
def f(x, y): return k * y

if __name__ == "__main__":
    '''
    This is a test to solve the simple system y' = ky
    with initial value (0, 1)
    '''


    solved = euler((0, 1))
    print(solved)

    lnsp = np.linspace(0, final_time, 1000)
    plt.plot(lnsp, solution(lnsp), c = (0., .1, .3), label = "Analytical solution")
    plt.xlabel("t")
    plt.yscale("log")

    plt.scatter(* zip(* solved), c = (.6, .2, .1), marker = "x", label = "Numerical solution")

    plt.legend()
    plt.show()