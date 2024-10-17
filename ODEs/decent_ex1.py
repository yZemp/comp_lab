import numpy as np
import matplotlib.pyplot as plt

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 100 # Default length of approximation


def unpack(arr):
    '''
    Receives array or array-like of n-dim vectors
    returns n arrays of coordinates
    '''
    coords = [[] for i in range(len(arr[0]))]

    for i in range(len(arr[0])):
        for j, el in enumerate(arr):
            # print(i, j)
            coords[i].append(arr[j][i])

    return coords


# Analytical solution
def solution(x):
    return np.sin(x)


def f(x, y, z, h):
    '''
    This represent phi of euler's method
    This represent known element of the differential equation

    NOT TO BE CONFUSED WITH _F(Y)

    NOTE: making argument like that instead of (Y) or (Y, h) is for simplicity's sake
    '''

    return - y

def rk2_f(x, y, z, h):
    '''
    This represent phi of rk2's method
    
    NOTE: making argument like that instead of (Y) or (Y, h) is for simplicity's sake
    '''

    return f(x + .5 * h, y + .5 * h * f(x, y, z, h), z, h)


def _F(Y, phi, h):
    '''
    This represent a vectorial function
    accepts a vector and returns a vector
    its the change rate of Y for every step

    NOT TO BE CONFUSED WITH f(x, y)
    '''

    return np.array([1, Y[2], phi(Y[0], Y[1], Y[2], h)], dtype = np.float128)


def _step(Y, h, phi):
    '''
    This is a "stepper" for one-step methods
    returns Y + h * _F(Y)
    where 
        _F is a generic function of Y that approximate the solution
            _F is a function of phi
            NOTE: with phi = f we have Euler's method
            NOTE: with phi = rk2_f we have RK2 method
        Y can be decomposed in (x, y, z)
            where:
                x = time
                y = actual y coordinate of the point
                z = derivative of y (only needed to calculate next Y)
    '''

    return(Y + h * _F(Y, phi, h))


def euler_method(Y0, h = H0, final_time = FINAL_TIME):
    '''
    Approximate the differential equation with known value f(x, y)
    Using Euler's method
    given:
        starting values (x0, y0, z0) = Y0 (this has to be an array-like)
        final_time = for how long to approx
        h = how little the steps (the smaller the more precise) 
    '''

    # This list will hold every step, starting with Y0
    steps = [np.array(Y0)]

    # Computing steps
    while steps[-1][0] < final_time:
        # NOTE: passing f to _F (through the stepper) means we are using Euler's method
        new_Y = _step(steps[-1], h, f)
        steps.append(new_Y)
    
    return steps


def rk2_method(Y0, h = H0, final_time = FINAL_TIME):
    '''
    Approximate the differential equation with known value f(x, y)
    Using RK2 algorithm
    given:
        starting values (x0, y0, z0) = Y0 (this has to be an array-like)
        final_time = for how long to approx
        h = how little the steps (the smaller the more precise) 
    '''

    # This list will hold every step, starting with Y0
    steps = [np.array(Y0)]

    # Computing steps
    while steps[-1][0] < final_time:
        # NOTE: passing rk2_f to _F (through the stepper) means we are using RK2's method
        new_Y = _step(steps[-1], h, rk2_f)
        steps.append(new_Y)
    
    return steps


if __name__ == "__main__":

    arrx = np.linspace(0., FINAL_TIME, 10_000)
    plt.plot(arrx, solution(arrx), c = (.1, .1, .1), label = "Analytical solution")

    ###############################################################################################
    # Euler method

    euler_solved = euler_method((0., 0., 1.))
    euler_coords = unpack(euler_solved)

    plt.plot(euler_coords[0], euler_coords[1], c = (.1, .7, .1), label = "Euler method")
    plt.legend()


    ###############################################################################################
    # RK2 method

    rk2_solved = rk2_method((0., 0., 1.))
    rk2_coords = unpack(rk2_solved)

    plt.plot(rk2_coords[0], rk2_coords[1], c = (.1, .1, .9), label = "RK2 method")
    plt.legend()
    plt.show()