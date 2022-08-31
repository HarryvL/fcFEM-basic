import numpy as np

from scipy.sparse import csc_matrix

row = np.array([0, 2, 2, 2, 0, 1, 2])

col = np.array([0, 2, 0, 1, 2, 2, 2])

data = np.array([1, -6, 2, 3, 4, 5, 6])

a = csc_matrix((data, (row, col)), shape=(3, 3))

for i in range(3):
    a[i,2] = 10
    a[2,i] = 10

print(a)
print()
print(a.toarray())