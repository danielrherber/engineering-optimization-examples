import numpy as np
from scipy.linalg import null_space
from numpy.random import rand

# data
A = np.array([[1,-1,0,0],
              [0,0,1,1]])
b = np.array([2,2])

# orthonormal basis for the null space of A
# Z = null(A)
# "rational" basis for the null space of A
Z = null_space(A)

# feasible point
xbar = np.array([3,1,0,2])

# check if it is feasible
print(A.dot(xbar) - b)

# add null space
v = rand(2,)*100
xnew = xbar + Z.dot(v)

# check that it is still feasible
print(A.dot(xnew)-b)