import numpy as np

# Solve system: Ux = b
# mat = U


# def inth_elem(mat, i, b, n):
#     if i == n: return b[i] / mat[i][i]
#     return x[i] = (b[i] - sum([mat[i][j] * x[j] for j in range(i + 1, n - 1)])) / mat[i][i]

def backward_sub(mat, b):
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
    return np.array([sum([mat[i][j] * vec[j] for j in range(len(vec))]) for i in range(len(mat))])





if __name__ == "__main__":
    mat = np.array([[2, 1, 1], [0, 1, -2], [0, 0, 1]])
    b = np.array([1, -1, 4])
    
    x = backward_sub(mat, b)
    print("Soluzione x = ", x)
    b_ihope = mat_vec_prod(mat, x)
    print("b (in teoria) = ", b)
