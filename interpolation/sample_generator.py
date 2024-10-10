import numpy as np
import scipy.stats as sts
import random as rand


if __name__ == "__main__":
    x = np.sort([100 * rand.random() for i in range(50)])
    y = [np.cos(el) + sts.norm.rvs(0, 5) for el in x]
    np.savetxt('x.txt', x)
    np.savetxt('y.txt', y)