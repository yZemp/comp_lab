import numpy as np
import matplotlib.pyplot as plt

######################################################################################
# Useful vars

H0 = .01 # Default precision of step
FINAL_TIME = 30 # Default length of approximation


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


def f(x, y):
    '''
    This represent known element of the differential equation

    NOT TO BE CONFUSED WITH _F(Y)
    '''
    return - y


def _F(Y):
    '''
    This represent a vectorial function
    accepts a vector and returns a vector
    its the change rate of Y for every step

    NOT TO BE CONFUSED WITH f(x, y)
    '''

    # NOTE: here we use f(x, y) as Phi

    return np.array([1, Y[2], f(Y[0], Y[1])], dtype = np.float128)


def _euler_step(Y, h):
    '''
    Euler "stepper"
    returns Y + h * f(x, y)
    where Y can be decomposed in (x, y, z)
    where:
        x = time
        y = actual y coordinate of the point
        z = y' derivative in y (only needed to calculate next Y)

    '''
    # print(Y + h * _F(Y), "\n")
    return(Y + h * _F(Y))


def euler_method(Y0, h = H0, final_time = FINAL_TIME):
    '''
    Approximate the differential equation with known value f(x, y)
    given:
        starting values (x0, y0, z0) = Y0 (this has to be an array-like)
        final_time = for how long to approx
        h = how little the steps (the smaller the more precise) 
    '''

    # This list will hold every step, starting with Y0
    steps = [np.array(Y0)]

    # Computing steps
    while steps[-1][0] < final_time:
        new_Y = _euler_step(steps[-1], h)
        steps.append(new_Y)
    
    return steps



if __name__ == "__main__":

    arrx = np.linspace(0., FINAL_TIME, 10_000)
    plt.plot(arrx, solution(arrx), label = "Analytical solution")

    ###############################################################################################
    # Euler method

    euler_solved = euler_method((0., 0., 1.))
    coords = unpack(euler_solved)

    plt.plot(coords[0], coords[1], label = "Euler method")
    plt.legend()
    plt.show()