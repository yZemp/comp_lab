import numpy as np
from matrix_utils import solve_linear_system, mat_vec_prod

if __name__ == "__main__":
    # Solve linear system
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = float)
    b = np.array([8, -2, 2])
    
    # NOTE: transforming the system in an equivalent one
    # ==> solution is the same
    x, mat, vec = solve_linear_system(mat, b)


    # print("Soluzione  = ", x)

    b = mat_vec_prod(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)