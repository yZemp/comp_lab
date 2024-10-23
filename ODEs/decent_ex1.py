import numpy as np
import matplotlib.pyplot as plt

######################################################################################
# Useful vars

H0 = .001 # Default precision of step
FINAL_TIME = 10 # Default length of approximation
START_VALS = (0., 0., 1.) # Starting values for x, y, z

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
    This represent known element of the differential equation

    incidentally:
    This represent phi of euler's method

    NOT TO BE CONFUSED WITH _F(Y)

    NOTE: passing (x, y, z, h) instead of (Y) or (Y, h) is for compatibility's sake
    '''

    return - y


def rk2_f(x, y, z, h):
    '''
    This represent phi of rk2's method
    
    NOTE: passing (x, y, z, h) instead of (Y) or (Y, h) is for compatibility's sake
    '''

    return f(x + .5 * h, y + .5 * h * f(x, y, z, h), z, h)


def rk4_f(x, y, z, h):
    '''
    This represent phi of rk2's method
    
    NOTE: passing (x, y, z, h) instead of (Y) or (Y, h) is for compatibility's sake
    '''

    k1 = f(x, y, z, h)
    k2 = f(x + .5 * h, y + .5 * h * k1, z, h)
    k3 = f(x + .5 * h, y + .5 * h * k2, z, h)
    k4 = f(x + h, y + h * k3, z, h)

    return (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def _step(Y, h, phi):
    '''
    This is a vectorial "stepper" for one-step methods
    returns Y + h * phi(Y)
    where 
        phi is a generic function of Y, h, f that approximate the solution
            NOTE: with phi = f we have Euler's method
            NOTE: with phi = rk2_f we have RK2 method
            NOTE: with phi = rk4_f we have RK4 method
        Y can be decomposed in (x, y, z)
            where:
                x = time
                y = actual y coordinate of the point
                z = derivative of y (only needed to calculate next Y)
    '''

    step = np.array([1, Y[2], phi(Y[0], Y[1], Y[2], h)], dtype = np.float128)

    return (Y + h * step)




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
        # NOTE: passing rk2_f to _F (through the stepper) means we are using the RK2 method
        new_Y = _step(steps[-1], h, rk2_f)
        steps.append(new_Y)
    
    return steps



def rk4_method(Y0, h = H0, final_time = FINAL_TIME):
    '''
    Approximate the differential equation with known value f(x, y)
    Using RK4 algorithm
    given:
        starting values (x0, y0, z0) = Y0 (this has to be an array-like)
        final_time = for how long to approx
        h = how little the steps (the smaller the more precise) 
    '''

    # This list will hold every step, starting with Y0
    steps = [np.array(Y0)]

    # Computing steps
    while steps[-1][0] < final_time:
        # NOTE: passing rk4_f to _F (through the stepper) means we are using the RK4 method
        new_Y = _step(steps[-1], h, rk4_f)
        steps.append(new_Y)
    
    return steps



def error_comparison(h0 = 1, hf = .01, time = FINAL_TIME):
    '''
    Plot errors as a function of h

    Different methods of computing errors:
        Error evaluated at final time (very random)
        Sum of all errors since t0 to final time 
    '''

    arrh = np.logspace(h0, hf, endpoint = True, base = 2, num = 1_000) - 1
    # arrh = np.linspace(h0, hf, endpoint = True, num = 1_000)
    # print(arrh)


    euler_err = []
    rk2_err = []
    rk4_err = []

    euler_err_fixed = []
    rk2_err_fixed = []
    rk4_err_fixed = []

    for h in arrh:
        euler_solved = euler_method(START_VALS, h = h)
        euler_coords = unpack(euler_solved)

        rk2_solved = rk2_method(START_VALS, h = h)
        rk2_coords = unpack(rk2_solved)

        rk4_solved = rk4_method(START_VALS, h = h)
        rk4_coords = unpack(rk4_solved)

        # TODO: find a better method of evaluation for errors
        
        euler_err.append(np.sum([abs(solution(euler_coords[0][i]) - euler_coords[1][i]) for i in range(len(euler_coords[0]))]))
        rk2_err.append(np.sum([abs(solution(rk2_coords[0][i]) - rk2_coords[1][i]) for i in range(len(rk2_coords[0]))]))
        rk4_err.append(np.sum([abs(solution(rk4_coords[0][i]) - rk4_coords[1][i]) for i in range(len(rk4_coords[0]))]))
        
        euler_err_fixed.append(abs(solution(euler_coords[0][-1]) - euler_coords[1][-1]))
        rk2_err_fixed.append(abs(solution(rk2_coords[0][-1]) - rk2_coords[1][-1]))
        rk4_err_fixed.append(abs(solution(rk4_coords[0][-1]) - rk4_coords[1][-1]))


    fig, ax = plt.subplots(2)
    
    ax[0].plot(arrh, euler_err, c = (.1, .7, .1), label = "Euler")
    ax[0].plot(arrh, rk2_err, c =  (.1, .1, .9), label = "RK2")
    ax[0].plot(arrh, rk4_err, c =  (.8, .1, .3), label = "RK4")

    ax[0].set_title("Error sum")
    # ax[0].set_yscale("log")
    # ax[0].set_xscale("log")
    ax[0].xaxis.set_inverted(True)

    ax[0].legend()



    ax[1].set_title("Error at fixed time")
    ax[1].plot(arrh, euler_err_fixed, c = (.1, .7, .1), label = "Euler")
    ax[1].plot(arrh, rk2_err_fixed, c =  (.1, .1, .9), label = "RK2")
    ax[1].plot(arrh, rk4_err_fixed, c =  (.8, .1, .3), label = "RK4")

    # ax[1].set_yscale("log")
    # ax[1].set_xscale("log")
    ax[1].xaxis.set_inverted(True)

    ax[1].legend()


    


    plt.show()




if __name__ == "__main__":

    arrx = np.linspace(0., FINAL_TIME, 10_000)
    plt.plot(arrx, solution(arrx), c = (.1, .1, .1), marker = "", label = "Analytical solution")


    ###############################################################################################
    # Euler method

    euler_solved = euler_method(START_VALS)
    euler_coords = unpack(euler_solved)

    plt.plot(euler_coords[0], euler_coords[1], marker = "", c = (.1, .7, .1), label = "Euler method")
    plt.legend()


    ###############################################################################################
    # RK2 method

    rk2_solved = rk2_method(START_VALS)
    rk2_coords = unpack(rk2_solved)

    plt.plot(rk2_coords[0], rk2_coords[1], marker = "", c = (.1, .1, .9), label = "RK2 method")
    plt.legend()


    ###############################################################################################
    # RK4 method

    rk4_solved = rk4_method(START_VALS)
    rk4_coords = unpack(rk4_solved)

    plt.plot(rk4_coords[0], rk4_coords[1], marker = "", c = (.8, .1, .3), label = "RK4 method")
    plt.legend()




    plt.show()


    ###############################################################################################
    # Errors

    error_comparison()