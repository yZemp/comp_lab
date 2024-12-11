import numpy as np
import matplotlib.pyplot as plt
from roots_lib import bisect, newton_raphson, secant, all_roots
from interpolation.Polynomial_classes import Polynomial

colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5), (.1, .2, .9)]


def f(x): return .5 + x - np.power(x, 2)
def g(x): return 1 - 2 * x

def f2(x): return np.power(x, 2)
def g2(x): return 2 * x

def f3(x): return np.power(x, 3) - 2 * x + 2
def g3(x): return 3 * np.power(x, 2) - 2

def f4(x): return np.power(x, 3) - 2 * np.power(x, 2) - 11 * x + 12
def g4(x): return 3 * np.power(x, 2) - 4 * x - 11

def f5(x): return np.sqrt(x + 1) * np.power(np.cos(x / 2), 3)
def g5(x): return np.power(np.cos(x / 2), 3) / (2 * np.sqrt(x + 1)) - (3 * np.power(np.cos(x / 2), 2) * np.sin(x / 2) * np.sqrt(x + 1)) / 2


def ex9_1():
    interval = [.8, 1.6]
    root_bis = bisect(*interval, f, .01)
    root_nr = newton_raphson(interval[0], f, g, .01)
    root_s = secant(*interval, f, .01)
    print("Bisection root:\t", root_bis)
    print("Newton-Raphson root:\t", root_nr)
    print("Secant root:\t", root_s)

    arrx = np.linspace(*interval, num = 10_000)
    fig, ax = plt.subplots()
    ax.plot(arrx, f(arrx), color = colors[0], label = "Function")
    ax.scatter(root_bis, f(root_bis), color = colors[1], marker = "x", label = "Root bisection")
    ax.scatter(root_nr, f(root_nr), color = colors[2], marker = "x", label = "Root bisection")
    ax.scatter(root_s, f(root_s), color = colors[3], marker = "x", label = "Root secant")
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')

    fig.set_size_inches((20, 10), forward = True)
    plt.legend()
    plt.savefig("roots_graphs/roots.png")
    plt.show()

def ex9_2():
    interval = [.8, 1.6]
    iternum = [i for i in range(1, 50)]
    root_bis = []
    root_nr = []
    root_s = []

    fig, ax = plt.subplots(2)

    for it in iternum:
        tmp_bis = bisect(*interval, f, max_iter = it)
        ax[0].scatter(tmp_bis, .1, marker = "x", color = (.4, .7 + .3 * it / 50, .1))
        root_bis.append(abs(f(tmp_bis)))

        tmp_nr = newton_raphson(interval[0], f, g, max_iter = it)
        ax[0].scatter(tmp_nr, .2, marker = "x", color = (.7 + .3 * it / 50, .1, .4))
        root_nr.append(abs(f(tmp_nr)))

        tmp_s = secant(*interval, f, max_iter = it)
        ax[0].scatter(tmp_s, .3, marker = "x", color = (.1, .4, .7 + .3 * it / 50))
        root_s.append(abs(f(tmp_s)))


    arrx = np.linspace(1.25, 1.45, num = 10_000)
    ax[0].plot(arrx, f(arrx), color = colors[0], label = "Function")
    ax[0].suptitle = "Roots convergence"

    ax[0].axhline(y=0, color='k')
    ax[0].legend()

    ax[1].plot(iternum, root_bis, color = colors[1], label = "Root convergence (bisection)")
    ax[1].plot(iternum, root_nr, color = colors[2], label = "Root convergence (Newton-Raphson)")
    ax[1].plot(iternum, root_s, color = colors[3], label = "Root convergence (secant)")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("Iteration")
    ax[1].set_ylabel("Error")

    fig.set_size_inches((20, 10), forward = True)
    ax[1].legend()
    plt.savefig("roots_graphs/convergence.png")
    plt.show()

def ex10_1():
    interval = [-.8, .8]
    iternum = [i for i in range(1, 50)]
    root_nr = []

    fig, ax = plt.subplots(2)

    for it in iternum:
        tmp_nr = newton_raphson(interval[0], f2, g2, max_iter = it)
        ax[0].scatter(tmp_nr, -.1, marker = "x", color = (.7 + .3 * it / 50, .1, .4))
        root_nr.append(abs(f2(tmp_nr)))


    arrx = np.linspace(*interval, num = 10_000)
    ax[0].plot(arrx, f2(arrx), color = colors[0], label = "Function")
    ax[0].suptitle = "Roots convergence"

    ax[0].axhline(y=0, color='k')
    ax[0].legend()

    ax[1].plot(iternum, root_nr, color = colors[2], label = "Root convergence (Newton_Raphson)")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("Iteration")
    ax[1].set_ylabel("Error")

    fig.set_size_inches((20, 10), forward = True)
    ax[1].legend()
    plt.savefig("roots_graphs/slow.png")
    plt.show()

def ex10_2(osc: bool):
    interval = [-2, 2]
    iternum = [i for i in range(1, 50)]
    root_nr = []
    root_s = []

    fig, ax = plt.subplots(2)

    for it in iternum:
        tmp_nr = newton_raphson(0, f3, g3, max_iter = it, fix_oscillations = osc)
        ax[0].scatter(tmp_nr, -.1, marker = "x", color = (.7 + .3 * it / 50, .1, .4))
        root_nr.append(abs(f3(tmp_nr)))
        
        tmp_s = secant(0, .001, f3, max_iter = it)
        ax[0].scatter(tmp_s, -.1, marker = "x", color = (.1, .4, .7 + .3 * it / 50))
        root_s.append(abs(f3(tmp_s)))


    arrx = np.linspace(*interval, num = 10_000)
    ax[0].plot(arrx, f3(arrx), color = colors[0], label = "Function")
    ax[0].suptitle = "Roots convergence"

    ax[0].axhline(y=0, color='k')
    ax[0].legend()

    ax[1].plot(iternum, root_nr, color = colors[2], label = "Root convergence (Newton_Raphson)")
    ax[1].plot(iternum, root_s, color = colors[3], label = "Root convergence (secant)")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("Iteration")
    ax[1].set_ylabel("Error")

    fig.set_size_inches((20, 10), forward = True)
    ax[1].legend()
    plt.savefig(f"roots_graphs/oscillating_{osc}.png")
    plt.show()

def ex10_3(starting_vals):
    interval = [-2, 2]
    iternum = [i for i in range(1, 50)]
    root_nr = [[] for s in starting_vals]

    fig, ax = plt.subplots(2)

    for i, s in enumerate(starting_vals):
        for it in iternum:
            tmp_nr = newton_raphson(s, f4, g4, max_iter = it)
            ax[0].scatter(tmp_nr, -.1, marker = "x", color = colors[i])
            root_nr[i].append(abs(f3(tmp_nr)))
        
        ax[1].plot(iternum, root_nr[i], color = colors[i], label = f"Newton_Raphson startval: {i + 1}")
            

    arrx = np.linspace(*interval, num = 10_000)
    ax[0].plot(arrx, f3(arrx), color = colors[0], label = "Function")
    ax[0].suptitle = "Roots convergence"

    ax[0].axhline(y=0, color='k')
    ax[0].legend()

    ax[1].set_yscale("log")
    ax[1].set_xlabel("Iteration")
    ax[1].set_ylabel("Error")

    fig.set_size_inches((20, 10), forward = True)
    ax[1].legend()
    plt.savefig(f"roots_graphs/sensitivity.png")
    plt.show()

def ex10_op():
    interval = [0, 2 * np.pi]
    iternum = [i for i in range(1, 50)]
    root_bis = []
    root_nr = []
    root_s = []

    fig, ax = plt.subplots(2)

    for it in iternum:
        tmp_bis = bisect(*interval, f5, max_iter = it)
        ax[0].scatter(tmp_bis, -.1 * 5, marker = "x", color = (.4, .7 + .3 * it / 50, .1))
        root_bis.append(abs(f5(tmp_bis)))

        tmp_nr = newton_raphson(1, f5, g5, max_iter = it)
        ax[0].scatter(tmp_nr, -.2 * 5, marker = "x", color = (.7 + .3 * it / 50, .1, .4))
        root_nr.append(abs(f5(tmp_nr)))

        tmp_s = secant(1, 1 + .001, f5, max_iter = it)
        ax[0].scatter(tmp_s, -.3 * 5, marker = "x", color = (.1, .4, .7 + .3 * it / 50))
        root_s.append(abs(f5(tmp_s)))


    arrx = np.linspace(0, 2 * np.pi, num = 10_000)
    ax[0].plot(arrx, f5(arrx), color = colors[0], label = "Function")
    ax[0].suptitle = "Roots convergence"

    ax[0].axhline(y=0, color='k')
    ax[0].legend()

    ax[1].plot(iternum, root_bis, color = colors[1], label = "Root convergence (bisection)")
    ax[1].plot(iternum, root_nr, color = colors[2], label = "Root convergence (Newton-Raphson)")
    ax[1].plot(iternum, root_s, color = colors[3], label = "Root convergence (secant)")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("Iteration")
    ax[1].set_ylabel("Error")

    fig.set_size_inches((20, 10), forward = True)
    ax[1].legend()
    plt.savefig("roots_graphs/optional.png")
    plt.show()

def ex11_leg():
    coeff = np.array([-63, 0, 3465, 0, -30030, 0, 90090, 0, -109395, 0, 46189]) / 256
    poly = Polynomial(coeff)
    eigvals, teigvals = all_roots(poly, N = 10_000)
    
    print(eigvals, "\n", teigvals)
    print("Close enough:\t", np.allclose(eigvals, teigvals))

def ex11_her():
    coeff = np.array([-120, 0, 720, 0, -480, 0, 64])
    poly = Polynomial(coeff)
    eigvals, teigvals = all_roots(poly, N = 10_000)
    
    print(eigvals, "\n", teigvals)
    print("Close enough:\t", np.allclose(eigvals, teigvals))


if __name__ == "__main__":
    starting_vals = [2.352837350, 2.352836327, 2.352836323]
    
    # ex9_1()
    # ex9_2()
    # ex10_1()
    # ex10_2(False)
    # ex10_2(True)
    # ex10_3(starting_vals)
    # ex10_op()
    ex11_leg()
    ex11_her()