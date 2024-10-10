import numpy as np
import matplotlib.pyplot as plt

def nth_newton_poly(j, n):
    # n is nth-chebyshev node
    if j == 0: return lambda x: 1

    def poly(x):
        return np.prod([(x - nth_chebyshev_node(i, n)) for i in range(0, j)])
     
    return poly

def nth_chebyshev_node(j, n):
    return - np.cos(j * np.pi / (n - 1))

# def nth_chebyshev_node(j, n):
#     def node(n):
#         return - np.cos(j * np.pi / (n - 1))
#     return node


if __name__ == "__main__":
    # Starting data
    x = [0, 10, 15, 20, 22.5, 30] # i=0,1,2,3,4,5
    f = [0, 227.04, 362.78, 517.35, 602.97, 901.67]

    lnsp = np.linspace(0, 30, 1_000)
    poly_4 = nth_newton_poly(2, 2)
    poly_values = np.array([poly_4(x) for x in lnsp])
    # print(poly_4(0))
    plt.plot(lnsp, poly_values, c = "#a909e3")

    plt.scatter(x, f, marker = "x", c = "#040404")



    plt.show()