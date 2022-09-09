import numpy as np
from scipy.optimize import NonlinearConstraint,BFGS,minimize
from matplotlib import pyplot as plt,patches

def wire_amount(x):

    '''
    Objective function
    '''

    # extract
    x0 = x[0]; y0 = x[1]
    x1 = x[2]; y1 = x[3]
    x2 = x[4]; y2 = x[5]
    x3 = x[6]; y3 = x[7]
    x4 = x[8]; y4 = x[9]

    # calculate individual wire lengths
    w1 = np.sqrt((x1-x0)**2 + (y1-y0)**2)
    w2 = np.sqrt((x2-x0)**2 + (y2-y0)**2)
    w3 = np.sqrt((x3-x0)**2 + (y3-y0)**2)
    w4 = np.sqrt((x4-x0)**2 + (y4-y0)**2)

    # total wire length
    w = w1 + w2 + w3 + w4

    return w


def building_constraints(x):

    '''
    Constraints
    '''

    # extract
    x0 = x[0]; y0 = x[1]
    x1 = x[2]; y1 = x[3]
    x2 = x[4]; y2 = x[5]
    x3 = x[6]; y3 = x[7]
    x4 = x[8]; y4 = x[9]

    # initialize
    c = np.zeros((10,))

    # building 1
    c[0] = (x1-1)**2 + (y1-4)**2 - 4

    # building 2
    c[1] = (x2-9)**2 + (y2-5)**2 - 1

    # building 3
    c[2] = 2 - x3
    c[3] = x3 - 4
    c[4] = -3 - y3
    c[5] = y3 + 1

    # building 4
    c[6] = 6 - x4
    c[7] = x4 - 8
    c[8] = -2 - y4
    c[9] = y4 - 2

    return c

# initial point
X0 = np.array([4,2,1,4,9,5,3,-2,7,0])

# upper and lower bound constraints for the nonlinear constraints to be of the form
# -inf <= c(x) <= 0

NL_lb = -np.ones((10,))*np.inf # define lower bound constraints
NL_ub = np.zeros((10,)) # define upper bound

# create the instance of the nonlinear constraints
nonlinear_constraint = NonlinearConstraint(building_constraints,NL_lb,NL_ub,jac = '2-point',hess = BFGS())

'''
optimization options

Check out https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
for a more comprehensive set of options.
Unfortunately not all of them can handle non-linear constraints.

'''
opt_method = 'COBYLA'
# opt_method = 'trust-constr'
# opt_method = 'SLSQP'

# solve the optimization problem
res = minimize(wire_amount,X0,method = opt_method,constraints = nonlinear_constraint,tol = 1e-8)

# extract results
x = res.x
F = res.fun

# extract
x0 = x[0]; y0 = x[1]
x1 = x[2]; y1 = x[3]
x2 = x[4]; y2 = x[5]
x3 = x[6]; y3 = x[7]
x4 = x[8]; y4 = x[9]

# initialize plot
fig = plt.figure()
ax = fig.add_subplot(111)

# plot buildings
circle1 = patches.Circle((1,4), radius=2, color='red') # building 1
circle2 = patches.Circle((9,5), radius=1, color='red') # building 2
rect1 = patches.Rectangle((2, -3), 2, 2, color='yellow') # building 3
rect2 = patches.Rectangle((6, -2), 2, 4, color='yellow') # building 4
ax.add_patch(rect1)
ax.add_patch(circle1)
ax.add_patch(rect2)
ax.add_patch(circle2)

ax.plot([x1,x0],[y1,y0],'-')
ax.plot([x2,x0],[y2,y0],'-')
ax.plot([x3,x0],[y3,y0],'-')
ax.plot([x4,x0],[y4,y0],'-')

plt.xlim([-1, 10])
plt.ylim([-3, 7])
plt.show()