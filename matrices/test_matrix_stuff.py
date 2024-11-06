import numpy as np
from matrix_utils import solve_linear_system, get_inverse, mat_vec_prod, LU_decompose, determinant

clear_space = "\n----------------------------------------------------------------\n"

if __name__ == "__main__":
    print(clear_space)

    # EX 1
    print("EX 1")

    # Solve linear system
    mat = np.array([[2, 1, 1], [0, 1, -2], [0, 0, 1]], dtype = float)
    b = np.array([1, -1, 4], dtype = float)
    
    x, _, _ = solve_linear_system(mat, b)

    print("Solution  = ", x)

    b = mat_vec_prod(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)

    print(clear_space)



    #EX 2
    print("EX 2")

    # Solve linear system
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = float)
    b = np.array([8, -2, 2], dtype = float)
    
    # NOTE: transforming the system in an equivalent one
    # ==> solution is the same
    x, mat, vec = solve_linear_system(mat, b)


    print("Solution  = ", x)

    b = mat_vec_prod(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)

    print(clear_space)



    #EX 2.1
    print("EX 2.1")

    # Invert matrix
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = float)
    
    inverse = get_inverse(mat)

    print("Inverse:\n", inverse, "\n")

    print("True inverse (numpy):\n", np.linalg.inv(mat), "\n")

    print(clear_space)



    #EX 2.2
    print("EX 2.2")

    # Solve linear system (needing partial pivoting)
    mat = np.array([[2, 1, 1], [2, 1, -4], [1, 2, 1]], dtype = float)
    b = np.array([8, -2, 2], dtype = float)
    
    x, mat, vec = solve_linear_system(mat, b)


    print("Solution  = ", x)

    b = mat_vec_prod(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)

    print(clear_space)




    #EX 2.3
    print("EX 2.3")

    # Solve linear system (needing partial pivoting)
    mat = np.array([[10e-20, -1], [1, 1]], dtype = np.float128)
    b = np.array([1, 2], dtype = np.float128)
    
    x, mat, vec = solve_linear_system(mat, b)


    print("Solution  = ", x)

    b = mat_vec_prod(mat, x)
    # b = np.dot(mat, x)
    print("Inverse operation (should yield b of the upper triangular system) = ", b)
    # FIXME: inverse operation is rekt

    print(clear_space)




    #EX 3
    print("EX 3")

    # Perform LU decomposition
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = float)
    
    L, U = LU_decompose(mat)


    print("L:\n", L)
    print("U:\n", U)

    check = np.dot(L, U)

    print("\nShould be original mat:")
    print(check)

    print("\nCalculating determinant using LU decomposition:\n")
    print(determinant(mat))

    print(clear_space)