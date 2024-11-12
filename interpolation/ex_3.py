import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
from matrices.matrix_utils import solve_linear_system
from ex_2 import chebyshev_nodes
from first_order_spline import linear_spline
from second_order_spline import quadratic_spline

from Polynomial_classes import Polynomial
from my_interplib import *

np.set_printoptions(linewidth=np.inf)

#######################################################################
# VARS

INTERVAL = (-.9, .9)
colors = [(.5, .1, .8), (.5, .8, .1), (.8, .1, .5)]
