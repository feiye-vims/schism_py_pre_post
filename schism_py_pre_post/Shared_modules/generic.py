def BinA(A=None, B=None):
    import numpy as np
    if A is None and B is None:
        A = np.array([3,5,7,1,9,8,6,6])
        B = np.array([3,1,5,8,6])
    else:
        A = np.array(A)
        B = np.array(B)

    index = np.argsort(A)
    sorted_A = A[index]
    sorted_index = np.searchsorted(sorted_A, B)

    Bindex = np.take(index, sorted_index, mode="raise")
    mask = A[Bindex] != B

    result = np.ma.array(Bindex, mask=mask)
    return result
