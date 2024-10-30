import numpy as np

def unpack(arr):
    '''
    Receives array or array-like of n-dim vectors
    returns n arrays of coordinates

    Example: 
        arr = [[x0, y0, z0], [x1, y1, z1], [x2, y2, z2]]
        returns: [[x0, x1, x2], [y0, y1, y2], [x2, y2, z2]]
    '''
    coords = [[] for i in range(len(arr[0]))]

    for i in range(len(arr[0])):
        for j, el in enumerate(arr):
            # print(i, j)
            coords[i].append(arr[j][i])

    # print(coords)
    return coords

def unpack_V2(arr):
    '''
    Receives array or array-like of vectors with time and [val1, val2]
    returns 2 arrays of time and y coords

    Example: 
        arr = [[t0, [y0, z0]], [t1, [y1, z1]], [t2, [y2, z2]]]
        returns: [[t0, t1, t2], [z0, z1, z2]]
    '''
    coords = [[] for i in range(len(arr[0][1]) + 1)]

    for j, el in enumerate(arr):
        # print(i, j)
        coords[0].append(arr[j][0])
        for k in range(len(arr[j][1])):
            coords[k + 1].append(arr[j][1][k])

    # print(coords)
    return coords