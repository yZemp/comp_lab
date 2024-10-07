import numpy as np


def to_upp_triangular(mat, b):
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
    '''

    for j, _ in enumerate(mat):
        for i in range(j + 1, len(mat)):
            print(f"\ni: {i}\tj: {j}\n")
            print(mat, b)
            # if mat[j][j] == 0: continue
            gauss_elim_coeff = mat[i][j] / mat[j][j]
            # TODO: can mat[j][j] ever be zero?
            print(mat[i], mat[j], gauss_elim_coeff)
            mat[i] = _gauss_elim(mat[i], mat[j], gauss_elim_coeff)
            b[i] = _gauss_elim(b[i], b[j], gauss_elim_coeff)

    print("Upper triangular form:")
    print(mat, b)

    return mat, b

def _gauss_elim(rowi, rowj, coeff):
    print("New row:\t ", rowi - coeff * rowj)
    return rowi - coeff * rowj


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
    print(x)
    for i in range(n - 1, -1, -1):
        print("i: ", i)
        # print("Somma\t", sum([j for j in range(i, i + 1)]))
        print("Somma\t", [j for j in range(i + 1, n - 1 + 2)])
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





if __name__ == "__main__":
    # Test
    mat = np.array([[2, 1, 1], [1, 1, -2], [1, 2, 1]], dtype = np.float128)
    b = np.array([8, -2, 2])

    mat, b = to_upp_triangular(mat, b)
    
    x = backward_sub(mat, b)
    print("Soluzione x = ", x)
    b = mat_vec_prod(mat, x)
    print("b (in teoria) = ", b)

