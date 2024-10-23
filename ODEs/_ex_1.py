import numpy as np
import matplotlib.pyplot as plt


final_time = 30


# Analytic solution
def solution(t): return np.sin(t)


# SOLVE METHODS

def euler(t0, theta0, theta_dot0, h = .001, final_time = final_time):
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


def RK4(t0, theta0, theta_dot0, h = .001, final_time = final_time):
    arr_t = [t0]
    arr_theta = [theta0]
    arr_theta_dot = [theta_dot0]

    domain = (t0, final_time)

    for i in range(int((domain[1] - domain[0]) / h)):
        new_t = arr_t[-1] + h


        k1 = f(arr_t[-1], arr_theta[-1])
        k2 = f(arr_t[-1] + .5 * h, arr_theta[-1] + .5 * h * k1)
        k3 = f(arr_t[-1] + .5 * h, arr_theta[-1] + .5 * h * k2)
        k4 = f(arr_t[-1] + h, arr_theta[-1] + h * k3)

        new_theta_dot = arr_theta_dot[-1] + h * (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        new_theta = arr_theta[-1] + h * arr_theta_dot[-1]

        arr_t.append(new_t)
        arr_theta.append(new_theta)
        arr_theta_dot.append(new_theta_dot)

    return [arr_t, arr_theta, arr_theta_dot]



def f(t, theta): return - theta


def main():
    '''
    Solving the system Theta_dot_dot = - Theta
    with initial value (Theta(0), Theta_dot(0)) = (0, 1)
    
    This is decomposed in:
        Theta_dot(t) = phi
        phi_dot(t) = - theta
    '''

    #########################################################################################
    # EULER METHOD
    euler_solved = euler(0., 0., 1.)
    # print(len(solved))

    fig, ax = plt.subplots(3)

    fig.set_size_inches(200, 100, forward = True)

    fig.suptitle("Euler's method", fontsize = 20)

    lnsp = np.linspace(0, final_time, 1000)
    ax[0].plot(lnsp, solution(lnsp), c = (0., .1, .3), label = "Correct $\Theta(t)$")
    plt.xlabel("t")
    # plt.yscale("log")

    xeuler = []
    yeuler = []
    for x, y, _ in zip(*euler_solved):
        xeuler.append(x)
        yeuler.append(y)

    ax[0].plot(xeuler, yeuler, c = (.1, .7, .3), label = "Euler h = .001")
    ax[0].legend()

    ax[1].plot(xeuler, solution(xeuler) - yeuler, c = (.9, .2, .5), label = f"Errors")
    ax[1].legend()


    # Perform multiple numeric integration with different h
    for h in [1, .8, .6, .4, .2, .1, .06, .04, .01, .005, .001]:
        euler_solved = euler(0., 0., 1., h = h)
        arrx = []
        arry = []
        for x, y, _ in zip(*euler_solved):
            arrx.append(x)
            arry.append(y)

        ax[2].plot(arrx, solution(arrx) - arry, c = (min([h, 1]), .2, .5), label = f"Errors h = {h}")
        ax[2].legend()

    ax[2].set_ylim([-50, 50])

    print("Euler computed")
    plt.show()
    fig.savefig("euler_harmonic_oscillator.png", dpi = fig.dpi, bbox_inches = "tight")
    plt.close()


    #########################################################################################
    # RK2 METHOD
    rk2_solved = RK2(0., 0., 1.)
    # print(len(solved))

    fig, ax = plt.subplots(3)
    fig.set_size_inches(200, 100, forward = True)
    fig.suptitle("RK2's method", fontsize = 20)

    lnsp = np.linspace(0, final_time, 1000)
    ax[0].plot(lnsp, solution(lnsp), c = (0., .1, .3), label = "Correct $\Theta(t)$")
    plt.xlabel("t")
    # plt.yscale("log")

    xrk2 = []
    yrk2 = []
    for x, y, _ in zip(*rk2_solved):
        xrk2.append(x)
        yrk2.append(y)

    ax[0].plot(xrk2, yrk2, c = (.1, .7, .3), label = "RK2 h = .001")
    ax[0].legend()

    ax[1].plot(xrk2, solution(xrk2) - yrk2, c = (.9, .2, .5), label = f"Errors")
    ax[1].legend()


    # Perform multiple numeric integration with different h
    for h in [1, .8, .6, .4, .2, .1, .06, .04, .01, .005, .001]:
        rk2_solved = RK2(0., 0., 1., h = h)
        arrx = []
        arry = []
        for x, y, _ in zip(*rk2_solved):
            arrx.append(x)
            arry.append(y)

        ax[2].plot(arrx, solution(arrx) - arry, c = (min([h, 1]), .2, .5), label = f"Errors h = {h}")
        ax[2].legend()

    ax[2].set_ylim([-50, 50])

    print("RK2 computed")
    plt.show()
    fig.savefig("rk2_harmonic_oscillator.png", dpi = fig.dpi, bbox_inches = "tight")
    plt.close()



    #########################################################################################
    # RK4 METHOD
    rk4_solved = RK4(0., 0., 1.)
    # print(len(solved))

    fig, ax = plt.subplots(3)
    fig.set_size_inches(200, 100, forward = True)
    fig.suptitle("RK4's method", fontsize = 20)

    lnsp = np.linspace(0, final_time, 1000)
    ax[0].plot(lnsp, solution(lnsp), c = (0., .1, .3), label = "Correct $\Theta(t)$")
    plt.xlabel("t")
    # plt.yscale("log")

    xrk4 = []
    yrk4 = []
    for x, y, _ in zip(*rk4_solved):
        xrk4.append(x)
        yrk4.append(y)

    ax[0].plot(xrk4, yrk4, c = (.1, .7, .3), label = "RK4 h = .001")
    ax[0].legend()

    ax[1].plot(xrk4, solution(xrk4) - yrk4, c = (.9, .2, .5), label = f"Errors")
    ax[1].legend()


    # Perform multiple numeric integration with different h
    for h in [1, .8, .6, .4, .2, .1, .06, .04, .01, .005, .001]:
        rk4_solved = RK4(0., 0., 1., h = h)
        arrx = []
        arry = []
        for x, y, _ in zip(*rk4_solved):
            arrx.append(x)
            arry.append(y)

        ax[2].plot(arrx, solution(arrx) - arry, c = (min([h, 1]), .2, .5), label = f"Errors h = {h}")
        ax[2].legend()

    ax[2].set_ylim([-50, 50])

    print("RK4 computed")
    plt.show()
    fig.savefig("rk4_harmonic_oscillator.png", dpi = fig.dpi, bbox_inches = "tight")
    plt.close()



    #########################################################################################
    # METHODS COMPARISON

    fig, ax = plt.subplots(2)
    fig.set_size_inches(200, 100, forward = True)
    fig.suptitle("Methods compared", fontsize = 20)

    ax[0].plot(xeuler, yeuler, label = "Euler")
    ax[0].plot(xrk2, yrk2, label = "RK2")
    ax[0].plot(xrk4, yrk4, label = "RK4")
    ax[0].legend()

    ax[1].plot(xeuler, yeuler - solution(xeuler), label = "Euler errors")
    ax[1].plot(xrk2, yrk2 - solution(xrk2), label = "RK2 errors")
    ax[1].plot(xrk4, yrk4 -  - solution(xrk4), label = "RK4 errors")
    ax[1].legend()
    

    plt.show()

if __name__ == "__main__":
    main()