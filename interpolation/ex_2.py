# from matrices.matrix_utils import get_inverse
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("~/Documents/Programming/comp_lab/")

#######################################################################
# VARS

INTERVAL = (-1, 1)




def runge(x):
    return 1 / (1 + 25. * x)



if __name__ == "__main__":
    arrx = np.linspace(*INTERVAL, num = 50)