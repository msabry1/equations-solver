import numpy as np


def lu_decomposition_crout(A):
    A = np.array(A, dtype=float)

    n, m = A.shape
    if n != m:
        raise ValueError("Input must be a non-empty 2D matrix")

    n = A.shape[0]

    # Check matrix is square
    if A.shape[1] != n:
        raise ValueError("Input matrix must be square")

    L = np.zeros(A, dtype=float)
    U = np.zeros(n, dtype=float) # identity matrix

    # we start by cols
    for j in range(n):
        U[j, j] = 1
        # Compute column j in L
        for i in range(j, n):
            L[i, j] = A[i, j] - np.sum(L[i, :j] * U[:j, j])

        if np.isclose(L[j, j], 0):
            raise ValueError("Matrix is singular and cannot be decomposed using LU decomposition")

        # Compute row j in U
        for i in range(j + 1, n):
            U[j, i] = (A[j, i] - np.sum(L[j, :j] * U[:j, i])) / L[j, j]

    return L, U


def verify_decomposition(A, L, U):
    AA = np.dot(L, U)

    # Check if reconstruction is close to original matrix
    return np.allclose(A, AA)


# Example usage
A = np.array([[4, -2, 1],
              [3, 6, -1],
              [2, -1, 5]])
L, U = lu_decomposition_crout(A)
is_valid = verify_decomposition(A, L, U)
# Expected output: True
print("Decomposition valid for A (3x3 matrix):", is_valid)


A = np.array([[2, -1],
              [-1, 2]])
L, U = lu_decomposition_crout(A)
is_valid = verify_decomposition(A, L, U)
# Expected output: True
print("Decomposition valid for A (2x2 matrix):", is_valid)


A = np.array([[1, 2],
              [3, 4],
              [5, 6]])
try:
    L, U = lu_decomposition_crout(A)
    is_valid = verify_decomposition(A, L, U)
except ValueError as e:
# Expected output: Error: Input matrix must be square
    print("Error:", e)


A = np.array([[0, 0],
              [0, 0]])
try:
    L, U = lu_decomposition_crout(A)
    is_valid = verify_decomposition(A, L, U)
except ValueError as e:
# Expected output: Error: Matrix is singular and cannot be decomposed using LU decomposition
    print("Error:", e)


A = np.array([[1, 2, 3],
              [2, 4, 6],
              [3, 6, 9]])
try:
    L, U = lu_decomposition_crout(A)
    is_valid = verify_decomposition(A, L, U)
except ValueError as e:
# Expected output: Error: Matrix is singular and cannot be decomposed using LU decomposition
    print("Error:", e)


A = np.eye(3)
L, U = lu_decomposition_crout(A)
is_valid = verify_decomposition(A, L, U)
# Expected output: True, L = I, U = I
print("Decomposition valid for Identity Matrix A:", is_valid)


A = np.random.rand(4, 4)
L, U = lu_decomposition_crout(A)
is_valid = verify_decomposition(A, L, U)
# Expected output: True
print("Decomposition valid for random 4x4 matrix:", is_valid)
