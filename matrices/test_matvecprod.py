import numpy as np
from backward_sub import backward_sub, mat_vec_prod

if __name__ == "__main__":
    # Solve 
    mat = np.array([[6, 563, 32, -27], [0, -14, 0, 14], [0, 0, 2, -1], [0, 0, 0, 87]])
    b = np.array([7, 2, 1, -5])
    
    x = backward_sub(mat, b)
    print("Soluzione x = ", x)
    b_ihope = mat_vec_prod(mat, x)
    print("b (in teoria) = ", b)
