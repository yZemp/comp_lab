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

    # partial pivoting to account for zero or close-to-zero elements in matjj:
    for j, _ in enumerate(mat):
        # print("-------------------------\nj cycle:\n", mat, "\n")
        print("Inizio:\n", mat, "\n")
        # Find row k > j with largest element in column
        index = 0
        max = 0
        for l in range(j, len(mat[j])):
            if mat[l][j] > max:
                max = mat[l][j]
                index = l
        if j == index: continue
        _gauss_swap(mat, j, index)
        _gauss_swap(b, j, index)
        print("Fine:\n", mat, "\n")

    # Gaussian elimination
    for j, _ in enumerate(mat):
        for i in range(j + 1, len(mat)):
            # print(f"\ni: {i}\tj: {j}\n")
            # print(mat, b)
            # if mat[j][j] == 0: continue
            gauss_elim_coeff = mat[i][j] / mat[j][j]
            
            mat[i] = _gauss_elim(mat[i], mat[j], gauss_elim_coeff)
            b[i] = _gauss_elim(b[i], b[j], gauss_elim_coeff)

    # print("Upper triangular form:")
    # print(mat, b)

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
    # print(id)

    inverse = np.array([[0. for j in range(len(mat))] for i in range(len(mat))], dtype = float)

    mat, id = _to_upp_triangular(mat, id)
    for i in range(len(mat)):
        # Solving system for every column of the identity
        col, _, _ = solve_linear_system(mat, [id_col[i] for id_col in id])
        for count, el in enumerate(col): inverse[count][i] = el

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
    
    mat, vec = _to_upp_triangular(mat, vec)
    # print("Equivalent upper triangular system:")
    # print(mat, vec)
    
    x = backward_sub(mat, vec)
    # print("Solution of the system:")
    # print(x)
    
    return x, mat, vec

if __name__ == "__main__":
    mat = np.array([[2, 1, 1], [2, 1, -4], [1, 2, 1]], dtype = float)
    b = np.array([8, -2, 2], dtype = float)
    x, _, _ = solve_linear_system(mat, b)
    print("X:")
    print(x) 

