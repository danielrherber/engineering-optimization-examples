import numpy as np
from scipy.optimize import linprog

# problem definition
f = np.array([100,150,200,400]) # profit coefficients
a1 = np.array([10,12,25,20]) # wood resourse coefficients
b1 = 5000 # amount of wood resource
a2 = np.array([2,4,8,12]) # labor resource coefficients
b2 = 1500 # amount of labor resource

LB = np.array([0,0,0,0]) # nonnegativity constraint
UB = np.array([None,None,None,None])

# combine bounds and reshape
bounds = np.vstack([LB,UB]); bounds = bounds.T

# combine A and B matrices
A = np.vstack([a1,a2]); B = np.vstack([b1,b2])

# specify method
method = 'highs'

# linprog options
options = {'disp':True}

# solve the optimization problem
res = linprog(-f,A_ub = A,b_ub = B,bounds = bounds,method = method, options = options)

# extract results
X = res.x # optimal solution
F = res.fun # optimal objective function value

# display optimal solution
print('')
print('bookshelves = {}'.format(X[0]))
print('')
print('cabinets with doors = {}'.format(X[1]))
print('')
print('tall cabinets with doors = {}'.format(X[2]))
print('')
print('fancy cabinets = {}'.format(X[3]))
print('')