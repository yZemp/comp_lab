import numpy as np


def _to_upp_triangular(mat, b):
    '''
    Accepts a system: Ux = b
    Where:
        U is a NxN upper triangular matrix
        x is unknown vector of length N
        b is known vector of length N

    It transform the generic NxN matrix
    to an upper triangular NxN matrix
    through Gaussian elimination
    and the b vector accordingly

    NOTE: Should also work if b is a Matrix
    '''

    # print(mat)
    # partial pivoting to account for zero or close-to-zero elements in matjj:
    for j, _ in enumerate(mat):
        # print("-------------------------\nj cycle:\n", mat, "\n")
        # print("Inizio:\n", mat, "\n")
        # Find row k > j with largest element in column
        index = j
        max_val = abs(mat[j][j])

        for l in range(j, len(mat[j])):
            if mat[l][j] > max_val:
                max_val = mat[l][j]
                index = l
        
        if j != index:    
            _gauss_swap(mat, j, index)
            _gauss_swap(b, j, index)
        
        # print("Fine:\n", mat, "\n")

    # Gaussian elimination
        for i in range(j + 1, len(mat)):
            # print(f"\ni: {i}\tj: {j}\n")
            # print(mat, b)
            # if mat[j][j] == 0: continue
            gauss_elim_coeff = mat[i][j] / mat[j][j]
            
            mat[i] = _gauss_elim(mat[i], mat[j], gauss_elim_coeff)
            b[i] = _gauss_elim(b[i], b[j], gauss_elim_coeff)

    # print("Upper triangular form:")
    # print("\n", mat, "\n", b)

    return mat, b

def _gauss_elim(rowi, rowj, coeff):
    '''
    Perform gaussian substitution of rowi with rowj times a specific coefficient
    '''
    # print("New row:\t ", rowi - coeff * rowj)
    return rowi - coeff * rowj

def _gauss_swap(arr, i, j):
    '''
    Swaps rows i, j of a given matrix / vector
    '''
    # FIXME: nothing works here
    
    # print("SWAPPING", i, j)
    temp = np.copy(arr[i])
    arr[i] = np.copy(arr[j])
    arr[j] = np.copy(temp)
    return

def get_inverse(mat):
    '''
    Returns inverse of the NxN matrix mat
    '''
    # Identity
    id = np.array([[1. if i == j else 0. for j in range(len(mat))] for i in range(len(mat))], dtype = float)
    inverse = np.zeros_like(mat, dtype = float)

    # mat, id = _to_upp_triangular(mat, id)
    # print(new_mat, "\n", new_id, "\n")
    
    for i in range(len(mat)):
        # Solving system for every column of the identity
        col, _, _ = solve_linear_system(np.copy(mat), np.copy(id[:, i]))
        inverse[:, i] = col
        # print("Inverse:\n", inverse, "\n")

    return inverse

def backward_sub(mat, b):
    '''
    Solve system: Ux = b
    Where:
        U is a NxN upper triangular matrix
        x is unknown vector of length N
        b is known vector of length N
    '''

    if not len(mat) == len(mat[0]) == len(b): return -1
    n = len(b) - 1

    x = np.zeros(n + 1)

    # FORMULA FOR i-NTH ELEMENT OF X
    # x[i] = (b[i] - sum([mat[i][j] * x[j] for j in range(i + 1, n - 1)])) / mat[i][i]


    x[n] = (b[n]) / mat[n][n]
    # print(x)
    for i in range(n - 1, -1, -1):
        # print("i: ", i)
        # print("Somma\t", sum([j for j in range(i, i + 1)]))
        # print("Somma\t", [j for j in range(i + 1, n - 1 + 2)])
        # plus two is for range() implementation #)
        x[i] = (b[i] - sum([mat[i][j] * x[j] for j in range(i + 1, n - 1 + 2)])) / mat[i][i]
        # print([mat[i][j] * x[j] for j in range(i + 1, n - 1)])
        # print(x[i])

    return x


def mat_vec_prod(mat, vec):
    '''
    Standard product Matrix x Vector
    With an NxN matrix
    and N length vector
    '''
    return np.array([sum([mat[i][j] * vec[j] for j in range(len(vec))]) for i in range(len(mat))])


def solve_linear_system(mat, vec):
    # print("Linear system to be solved:")
    # print(mat, vec)
    
    new_mat, new_vec = _to_upp_triangular(np.copy(mat), np.copy(vec))
    # print("Equivalent upper triangular system:")
    # print(mat, vec)
    
    x = backward_sub(new_mat, new_vec)
    # print("Solution of the system:")
    # print(x)
    
    return x, new_mat, new_vec

def LU_decompose(mat):
    '''
    Decompose a NxN matrix A in two matrices L, U
    where:
        L is lower triangular with unitary diagonal
        U is upper triangular 
    '''
    
    L = np.array([[1. if j == i else 0. for j in range(len(mat))] for i in range(len(mat))], dtype = float)

    # Gaussian elimination
    for j, _ in enumerate(mat):
        for i in range(j + 1, len(mat)):
            # print(j, i)
            # print(f"\ni: {i}\tj: {j}\n")
            # print(mat, b)
            # if mat[j][j] == 0: continue
            gauss_elim_coeff = mat[i][j] / mat[j][j]
            
            mat[i] = _gauss_elim(mat[i], mat[j], gauss_elim_coeff)

            L[i][j] = gauss_elim_coeff

    return L, mat

def determinant(mat):
    '''
    Calculate determinant of NxN matrix using the LU decomposition method
    '''

    L, U = LU_decompose(mat)
    # Det(L) = 0
    # U is upper triangular
    return np.prod([row[i] for i, row in enumerate(U)])

if __name__ == "__main__":
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = float)
    print(determinant(mat))