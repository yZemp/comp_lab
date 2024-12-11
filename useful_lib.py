import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def _cround(z, threshold = 10e-8):
    if abs(z.imag) < threshold:
        return z.real
    return z
