import numpy as np
import matplotlib.pyplot as plt

def euler(init_value: tuple, h = .01):
    solved = [init_value]

    for i in range(1000):
        new_x = solved[-1][0] + h
        new_y = y(solved[-1][0]) + h * k *y(solved[-1][0])
        solved.append((new_x, new_y))

    return solved

def y(t): return np.exp(t)


k = 14

if __name__ == "__main__":
    '''
    This is a test to solve the simple system y' = ky
    with initial value (0, 1)
    '''


    solved = euler((0, 1))
    print(solved)

    lnsp = np.linspace(0, 10, 1000)
    plt.plot(lnsp, y(lnsp), c = (0., .1, .3), label = "Correct y(t)")
    plt.xlabel("t")
    # plt.yscale("log")

    plt.scatter(* zip(* solved), c = (.6, .2, .1), marker = "x", label = "Numerical solution")

    plt.legend()
    plt.show()