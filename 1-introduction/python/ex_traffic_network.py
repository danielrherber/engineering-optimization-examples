import numpy as np
from scipy.optimize import minimize ,LinearConstraint,Bounds


def travel_time(x,*p):

    # extract
    t = p[0]; a = p[1]; c = p[2]

    # calculate total travel time
    f = x*(t + a*(x/(np.ones((5,))-(x/c))))

    f = np.sum(f)

    return f

def calc_jac(x,*p):

    # extract
    t = p[0]; a = p[1]; c = p[2]

    t0 = x*a

    t1 = 1-x/c

    jac = t + 2*t0/t1 + t0*x/(t1*t1)/c

    return jac

# problem data
t = np.array([0.1,0.14,0.13,0.25,0.12]) # constant travel time [hour]
a = np.array([0.0001,0.0001,0.0002,0.0001,0.0001]) # travel rate density constant [hour^2/car]
c = np.array([2200,2200,1500,1000,2000]) # road capacities [car/hour]
X = 2000 # volume of incoming/outgoing cars [car/hour]
ep = 0.1 # no negative rates [car/hour]

# argument tuple for minimizer
p = tuple([t,a,c])

# initial guess (happens to be feasible)
X0 = np.array([X,0,0,X,0])

# bounds
LB = np.zeros((5,))
UB = c-ep

# equality constraints
Aeq = np.array([[1,1,0,0,0],
                [-1,0,1,1,0],
                [0,-1,-1,0,1],
                [0,0,0,1,1]])

beq = np.array([X,0,0,X])

'''
Check out https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.LinearConstraint.html#scipy.optimize.LinearConstraint
for more details on the LinearConstraint object
'''

lincon = LinearConstraint(Aeq,beq,beq)
bounds = Bounds(LB,UB)

# optimization method
opt_method = 'SLSQP'

# solve the optimization
res = minimize(travel_time,X0,args = p,method = opt_method,constraints = (lincon,),tol = 1e-8 )
