import numpy as np

N = 3
num = 1000

B = np.zeros((num, 2 * N))
D = np.zeros((num, N))
V = np.zeros((num, N * N))


for i in range(num):
    # 20 October 2021
    # https://stackoverflow.com/questions/10806790/generating-symmetric-matrices-in-numpy/27331415
    b = np.random.uniform(-100, 100, size=(N, N))
    b = (b + b.T)/2

    # eigenvalues are in ascending order, we need descending order
    dd, vv = np.linalg.eigh(b)

    for j in range(3):
        D[i, j] = dd[2-j]
        V[i, 3*j:3*(j+1)] = vv[:, 2-j]

    B[i, 0:3] = b[0, :]
    B[i, 3:5] = b[1, 1:]
    B[i, 5] = b[2, 2]

np.savetxt('B.asc', B, header=str(num), comments='')
np.savetxt('D.asc', D, header=str(num), comments='')
np.savetxt('V.asc', V, header=str(num), comments='')
