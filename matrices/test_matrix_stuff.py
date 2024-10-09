import numpy as np
from matrix_utils import solve_linear_system, get_inverse, mat_vec_prod

clear_space = "\n----------------------------------------------------------------\n"

if __name__ == "__main__":
    print(clear_space)

    # EX 1

    # Solve linear system
    mat = np.array([[2, 1, 1], [0, 1, -2], [0, 0, 1]], dtype = float)
    b = np.array([1, -1, 4])
    
    x, _, _ = solve_linear_system(mat, b)

    print("Solution  = ", x)

    b = mat_vec_prod(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)

    print(clear_space)



    #EX 2

    # Solve linear system
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = float)
    b = np.array([8, -2, 2])
    
    # NOTE: transforming the system in an equivalent one
    # ==> solution is the same
    x, mat, vec = solve_linear_system(mat, b)


    print("Solution  = ", x)

    b = mat_vec_prod(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)

    print(clear_space)



    #EX 3

    # Invert matrix
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = float)
    
    inverse = get_inverse(mat)

    print("Inverse:\n", inverse)

    print(clear_space)



    #EX 4

    # Solve linear system (needing partial pivoting)
    mat = np.array([[2, 1, 1], [2, 1, -4], [1, 2, 1]], dtype = float)
    b = np.array([8, -2, 2])
    
    x, mat, vec = solve_linear_system(mat, b)


    print("Solution  = ", x)

    b = mat_vec_prod(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)

    print(clear_space)