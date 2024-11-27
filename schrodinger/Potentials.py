import numpy as np
import matplotlib.pyplot as plt
from typing import Union

FloatOrArray = Union[float, np.ndarray]

class Potential():
    '''
    Different kind of potential objects
    '''

    def __init__(self, actual_func):
        self.actual_func = actual_func

    def __call__(self, x: FloatOrArray) -> float:
        return self.actual_func(x)

    def __str__(self):
        return f"Potential object:\n{self.actual_func}"


class Potential_simple(Potential):
    '''
    Simple is meant in the mathmatical sense of a simple function,
    defined by constants over descrete integrals.
    This kwargs is thought to be V0, V1, V2, ...
    and args is thought to be [a, b, c, d, ...] the interval of definition (else is 0)
    '''
    def __init__(self, args, **kwargs):
        def func(x):
            # If x is array
            if isinstance(x, np.ndarray):
                result = np.zeros_like(x)
                for i, (a, b) in enumerate(zip(args[:-1], args[1:])):
                    mask = (a <= x) & (x < b)  # Boolean mask
                    result[mask] = kwargs.get(f"V{i}", 0)
                return result
        
            # If x is a value
            for i, (a, b) in enumerate(zip(args[:-1], args[1:])):
                if a <= x < b:
                    return kwargs.get(f"V{i}", 0)

            return 0.
        
        super().__init__(func)


if __name__ == "__main__":
    potential = Potential_simple([-2., 1., 4., 6.], V0 = -10, V1 = 0, V2 = 20)
    # potential = Potential_simple([-2., 3.], V0 = -10)
    print(potential)
    arrx = np.linspace(-10, 10, num = 1000)
    plt.plot(arrx, potential(arrx), color = (.1, .1, .1), label = "Simple potential")
    plt.legend()
    plt.show()