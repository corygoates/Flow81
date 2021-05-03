import numpy as np

N = 2000
matrix = np.zeros((N,N))

for i in range(N):
    for j in range(N):
        matrix[i,j] = i*j
        matrix[i,j] = matrix[i,j]-i+j

for i in range(N):
    for j in range(N):
        print(matrix[i,j])