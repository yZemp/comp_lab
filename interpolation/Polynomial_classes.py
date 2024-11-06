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
    

    def __repr__(self):
        return f"P(x) = " + " + ".join([f"{elem} * x^({i})" for i, elem in enumerate(self.coefficients)])


    def evaluate(self, x):
        '''
        Evaluates the polynomial at a given point x.
        
        x: float or array-like
        return: polynomial value at x
        '''

        return np.sum([el * np.power(x, i) for i, el in enumerate(self.coefficients)], dtype = np.float128)



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
    
    def __repr__(self):
        return f"n_{self.order} = 1" + "".join([f"(x - {elem})" for elem in self.elements])


    def evaluate(self, x):
        '''
        Evaluates the Newton polynomial at a given point x.
        
        x: float or array-like
        return: Newton polynomial value at x
        '''

        if self.order == 0:
            return 1
        else:
            return np.prod([x - el for el in self.elements])




class Newton_interpolator:
    def __init__(self, dd, newton_poly):
        '''
        array: array-like

        Sets up a Newton polynomial with
            dd, newton_poly array
            dd are the devided differnces
            newton_poly are known points
        
            Polynomial is: P(x) = f_0 + f_01 (x - x0) + f_012 (x - x0) (x - x1) + f_0123 (x - x0) (x - x1) (x - x2) + ...
        '''

        try: dd = np.array(dd, dtype = np.float128)
        except: raise ValueError("Param should be an array-like (dd)")
        else:
            self.dd = dd
    
        try: newton_poly = np.array(newton_poly)
        except: raise ValueError("Param should be an array-like (newton_poly)")
        else:
            self.newton_poly = newton_poly
            # self.order = len(newton_poly) - 1


    def evaluate(self, x):
        '''
        Evaluates the interpolator at a given point x.
        
        x: float or array-like
        return: Interpolated value at x
        '''
        
        # print(self.newton_poly)
        return  np.sum([f * n.evaluate(x) for f, n in zip(self.dd, self.newton_poly)])



    
if __name__ == "__main__":
    pass