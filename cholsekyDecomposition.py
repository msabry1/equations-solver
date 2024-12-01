import numpy as np


def cholesky_decomposition(A):
    A = np.array(A, dtype=float)

    n, m = A.shape
    if n != m:
        raise ValueError("Matrix must be square")

    if not np.allclose(A, A.T):
        raise ValueError("Matrix must be symmetric")

    if not np.all(np.linalg.eigvals(A) >= 0):
        raise ValueError("Matrix must be positive definite")

    L = np.zeros_like(A)

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                sum = np.dot(L[i, :i], L[i, :i]) # sum of pow(L[i, j], 2)
                L[i, i] = np.sqrt(A[i, i] - sum)
            else:
                sum = np.dot(L[i, :j], L[j, :j]) # sum of L[i, j] * L[j, j]
                L[i, j] = (A[i, j] - sum) / L[j, j]


    return L, L.T

def verify_cholesky_decomposition(A, L, U):
    AA = L @ U
    return np.allclose(A, AA)


A = [[4, 12, -16],
     [12, 37, -43],
     [-16, -43, 98]]

L, U = cholesky_decomposition(A)

# Verify the decomposition
is_valid = verify_cholesky_decomposition(A, L, U)
print(L)
print(U)
print("Decomposition valid for A:", is_valid)
# Should output: True

A = [[25, 15, -5],
     [15, 18, 0],
     [-5, 0, 11]]

L, U = cholesky_decomposition(A)

# Verify the decomposition
is_valid = verify_cholesky_decomposition(A, L, U)
print("Decomposition valid for A:", is_valid)
# Should output: True

# Generate a random symmetric positive definite matrix
np.random.seed(42)
B = np.random.rand(3, 3)
A = B.T @ B  # This ensures A is symmetric and positive definite

L, U = cholesky_decomposition(A)

# Verify the decomposition
is_valid = verify_cholesky_decomposition(A, L, U)
print("Decomposition valid for random A:", is_valid)
# Should output: True

A = [[1, 2],
     [3, 4]]
try:
    L, U = cholesky_decomposition(A)
    is_valid = verify_cholesky_decomposition(A, L, U)
except ValueError as e:
    print("Error:", e)
# Should output: Error: Input matrix must be symmetric

A = [[1, 2],
     [2, 1]]
try:
    L, U = cholesky_decomposition(A)
    is_valid = verify_cholesky_decomposition(A, L, U)
except ValueError as e:
    print("Error:", e)
# Should output: Error: Input matrix is not positive definite
