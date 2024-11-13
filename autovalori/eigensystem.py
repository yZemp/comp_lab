import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from matrices.matrix_utils import mat_vec_prod

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

N_MAX = 30
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]


def _rayleigh_quotient(A, yn):
    num = np.vdot(np.transpose(yn), np.dot(A, yn))
    den = np.linalg.norm(yn)
    return num / den

def power_method(mat, N = N_MAX):
    '''
    Finds and returns the biggest eigenvalue of a square matrix
    '''

    x0 = np.array([row[0] for row in mat])
    y_n = [x0]
    for _ in range(N):
        new_y = np.dot(mat, y_n[-1])
        new_y_normalized = new_y / np.linalg.norm(new_y)
        y_n.append(new_y_normalized)
    
    approx = _rayleigh_quotient(mat, y_n[-1])
    return approx

if __name__ == "__main__":
    pass