import numpy as np

class Polynomial:
    def __init__(self, param):
        '''
        param: int or array-like

        If param is an int, sets up a polynomial of that order == param where all coefficients are 1
        If param is an array, sets up a polynomial with coefficients of the param array
        '''

        if isinstance(param, int):
            self.order = param
            self.coefficients = np.ones(param, dtype = np.float128)
        
        else:
            try: coeff = np.array(param, dtype = np.float128)
            except: raise ValueError("Param should be either an int (order) or an array-like (coefficients)")
            else:
                self.coefficients = coeff
                self.order = len(coeff)
    

    def evaluate(self, x):
        '''
        Evaluates the polynomial at a given point x.
        
        x: float or array-like
        return: polynomial value at x
        '''

        if self.order == 0:
            return 1
        return np.sum([el * np.power(x, i) for i, el in enumerate(self.coefficients)], dtype=np.float128)



class Newton_polynomial:
    def __init__(self, array):
        '''
        array: array-like

        Sets up a Newton polynomial with x_j elements of array
        '''
        try: elements = np.array(array, dtype = np.float128)
        except: raise ValueError("Param should be an array-like (x_j)")
        else:
            self.elements = elements
            self.order = len(elements)
    

    def evaluate(self, x):
        '''
        Evaluates the Newton polynomial at a given point x.
        
        x: float or array-like
        return: Newton polynomial value at x
        '''

        if self.order == 0:
            return 1
        return np.prod([x - el for el in self.elements])
