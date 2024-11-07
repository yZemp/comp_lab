import numpy as np


def chebyshev_nodes(N = 50):
    return [- np.cos((j * np.pi) / (N - 1)) for j in range(N)]


def runge(x):
    return 1 / (1 + 25. * np.power(x, 2))
