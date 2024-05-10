# Metode Numerik C

# Metode dekomposisi LU Gauss
def luGauss(A, b):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]
    for i in range(n):
        L[i][i] = 1.0
        for j in range(i, n):
            U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        for j in range(i+1, n):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    y = [0.0] * n
    x = [0.0] * n
    for i in range(n):
        y[i] = b[i] - sum(L[i][k] * y[k] for k in range(i))
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(U[i][k] * x[k] for k in range(i+1, n))) / U[i][i]

    return x

# Metode dekomposisi Crout
def crout(A, b):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1.0
        for j in range(i, n):
            U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        for j in range(i+1, n):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    y = [0.0] * n
    x = [0.0] * n
    for i in range(n):
        y[i] = b[i] - sum(L[i][k] * y[k] for k in range(i))
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(U[i][k] * x[k] for k in range(i+1, n))) / U[i][i]

    return x

# Metode dekomposisi Matriks Balikan
def matriksBalikan(A, b):
    A_inv = [[0.0] * len(A) for _ in range(len(A))]
    det = determinan(A)
    if det == 0:
        return "Matriks tidak memiliki invers karena determinan = 0"
    else:
        for i in range(len(A)):
            for j in range(len(A)):
                minor = [row[:j] + row[j+1:] for row in (A[:i]+A[i+1:])]
                A_inv[i][j] = (-1)**(i+j) * determinan(minor) / det

        x = [0.0] * len(A)
        for i in range(len(A)):
            x[i] = sum(A_inv[i][j] * b[j] for j in range(len(A)))
        return x

# Fungsi untuk menghitung determinan matriks
def determinan(A):
    if len(A) == 1:
        return A[0][0]
    elif len(A) == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    else:
        det = 0
        for c in range(len(A)):
            minor = [row[:c] + row[c+1:] for row in A[1:]]
            det += ((-1)**c) * A[0][c] * determinan(minor)
        return det

# Testing
J = [[4, 5, 6],
     [1, 3, 8],
     [3, 1, 2]]
K = [1, 2, 2]

print("Nama:", "Salman Arya Sandytia")
print("NIM:", 21120122140146)

print("\nMetode Dekomposisi LU Gauss:")
print(luGauss(J, K))
print("\nMetode Dekomposisi Crout:")
print(crout(J, K))
print("\nMetode Matriks Balikan:")
print(matriksBalikan(J, K))
