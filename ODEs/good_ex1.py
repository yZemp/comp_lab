import numpy as np
import matplotlib.pyplot as plt

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 20 # Default length of approximation
STEPS_NUMBER = 5_000 # Default steps number
START_VALS = (0., 0., 1.) # Starting values for t, z1, z2, ... , zN
gamma = .6
A = 1.4

def unpack(arr):
    '''
    Receives array or array-like of n-dim vectors
    returns n arrays of coordinates

    Example: 
        arr = [[x0, y0, z0], [x1, y1, z1], [x2, y2, z2]]
        returns: [[x0, x1, x2], [y0, y1, y2], [x2, y2, z2]]
    '''
    coords = [[] for i in range(len(arr[0]))]

    for i in range(len(arr[0])):
        for j, el in enumerate(arr):
            # print(i, j)
            coords[i].append(arr[j][i])

    # print(coords)
    return coords


# Analytical solution
def solution(t): return np.sin(t)



# Known term of ODE (in terms of t, z1, ... , zN)
def f(t, *args):
    '''
    Known term of the ODE (in terms of t, z1, ... , zN)
    args is expected to be z1, ... , zN
    '''

    return - args[0]
    # return - np.sin(args[0]) - gamma * args[1] + A * np.sin((2 / 3) * t)


def euler_phi(Y, h, known):
    return known(*Y)


def rk2_phi(Y, h, known):

    k1 = known(*Y)
    tmp = np.copy(Y)
    tmp[0] += .5 * h
    tmp[2] += .5 * h * k1
    tmp[1] += .5 * h * Y[2]

    return known(*tmp)



def _step(Y, h, phi, f):
    '''
    This is a vectorial "stepper" for one-step methods
    returns Y + h * phi(Y)
    where 
        phi is a generic function of Y, h, f that approximate the solution
        Y can be decomposed in (t, z1, ... , zN, zN_dot)
            where:
                t = time
                zN_d = actual y coordinate of the point
                zi (from 1 to N) = elements of first order system (only needed to calculate next Y)
    '''

    # Creating what to add every step 
    tmp = np.array(Y)[2:]
    step = np.concatenate((np.array([1.]), np.array([phi(Y, h = h, known = f) for el in tmp])))
    # FIXME: I don't really know how to fix the line above. Fuck.

    return (Y + h * step)




def solve(Y0, final_time = FINAL_TIME, h = H0, steps_number = -1, phi = euler_phi, f = f):
    '''
    Y0 should be an array-like of the form: [t_0, z1_0, ... , zN_0, zN_dot_0]

    Using:
        phi = f means Euler's method
    '''

    steps_number = int(final_time / h) if steps_number == -1 else steps_number

    steps = [Y0]

    for i in range(steps_number):
        new_step = _step(steps[-1], h, phi = phi, f = f)
        steps.append(new_step)

    return steps



if __name__ == "__main__":
    '''
    Solve a generic system:
        y_(N+1)dot = f(x, y, y_dot, ... , y_(N)dot)

    That can be rewritten in
        y_dot = z1
        z1_dot = z2
        ...
        zN_dot = f(t, z1, z2, ... , zN)

    For every step, Y is the vector (t, y_dot, z1_dot, z2_dot, ... , zN_dot) with
        t = time
        y_dot = height of point
        z_i_dot = first derivative of z_i for every i, or i-nth derivative of y
    '''

    arrx = np.linspace(0., STEPS_NUMBER * H0, 10_000)
    plt.plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution") # I DONT HAVE IT :(

    euler_values = solve(START_VALS, phi = euler_phi)
    euler_coords = unpack(euler_values)


    rk2_values = solve(START_VALS, phi = rk2_phi)
    rk2_coords = unpack(rk2_values)


    plt.plot(euler_coords[0], euler_coords[1], label = "Euler solution")
    plt.plot(rk2_coords[0], rk2_coords[1], label = "RK2 solution")
    plt.legend()
    plt.show()