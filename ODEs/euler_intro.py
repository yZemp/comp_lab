import numpy as np
import matplotlib.pyplot as plt

def euler(init_value: tuple, h = .01, domain = (0, 10)):
    solved = [init_value]

    for i in range(int((domain[1] - domain[0]) / h)):
        new_x = solved[-1][0] + h
        new_y = solved[-1][0] + h * f(solved[-1][0], solved[-1][1])
        solved.append((new_x, new_y))

    return solved


# Analytic solution

def solution(t): return np.exp(t)


k = 14
def f(x, y): return k * y

if __name__ == "__main__":
    '''
    This is a test to solve the simple system y' = ky
    with initial value (0, 1)
    '''


    solved = euler((0, 1))
    print(solved)

    lnsp = np.linspace(0, 10, 1000)
    plt.plot(lnsp, solution(lnsp), c = (0., .1, .3), label = "Analytical solution")
    plt.xlabel("t")
    plt.yscale("log")

    plt.scatter(* zip(* solved), c = (.6, .2, .1), marker = "x", label = "Numerical solution")

    plt.legend()
    plt.show()