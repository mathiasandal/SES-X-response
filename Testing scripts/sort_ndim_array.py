import numpy as np

a = np.ones([2, 2])


def my_func(mat):
    mat[0, 0] = 10

my_func(a)

print(a)



