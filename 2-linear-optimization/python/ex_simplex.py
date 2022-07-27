import numpy as np 
from scipy.optimize import linprog 
from scipy.linalg import solve


# simplex iteration
# NOTE: this assumes that you start at a basic feasible solution
def simplex_iteration(A,b,c,x):
    
    # display current value of x
    print('x = {}'.format(x))
    
    # initial values
    isUnbounded = False; isOptimal = False
    
    # get the number of variables, basic variables, and nonbasic variables
    nb,n = np.shape(A)
    nn = n-nb
    
    # determine the indices of nonbasic variables
    ind = list(np.argsort(x)[:nn])
    ind_ = np.arange(n)
    
    # create In and Ib as truth arrays
    In = np.array(n*[False])
    In[ind] = True
    Ib = ~In
    
    #Ib = list(Ib); In = list(In)
    
    # extract N and B matrices
    N = A[:,In]
    B = A[:,~In]
    
    # extract cN and cB
    cN = c[In]; cB = c[~In]
    
    # compute simplex multipliers
    y = solve(B.T,cB)
    
    # compute reduced cost vector
    chat = cN - np.dot(N.T,y)
    
    if all(chat >=0):
        print('Optimal!!!')
        isOptimal = True
        return x,isOptimal,isUnbounded
    
    # select entering variable        
    It_ = np.argmin(chat)
    It = ind_[In][It_]
    
    # step computation
    At = solve(B,A[:,It])
    bt = solve(B,b)
    
    # compute ratios
    Alpha = bt/At
    
    Alpha[Alpha<0] = np.inf
    
    # check if problem is unbounded
    if all(Alpha == np.inf):
        print('Unbounded!')
        isUnbounded = True 
        return x,isOptimal,isUnbounded
    
    # find smallest positive number for value of the entering variable 
    xt = np.min(Alpha)
    
    # update
    Z_ = np.vstack([np.eye(nn),solve(-B,N)])
    row_ = np.hstack([ind_[In],ind_[Ib]])
    
    # restore Z back to original ordering
    Z = np.zeros(np.shape(Z_))
    Z[row_,:] = Z_
    
    # step using correct column of Z
    x = x + xt*Z[:,It_]

    return x,isOptimal,isUnbounded
    
    

# test number (see below)
test = 4

if test == 1: # session 5 Class example
    A = np.array([[-2,1,1,0,0],
        [-1,2,0,1,0],
        [1,0,0,0,1]
        ])
    b = np.array([2,7,3])
    c = np.array([-1,-2,0,0,0])
    x = np.array([0,0,2,7,3])
    
elif test == 2: # Example for the Two-phase Method Continued (Phase 1)
    A = np.array([[3,2,0,0,1,0],
                  [2,-4,-1,0,0,1],
                  [4,3,0,1,0,0]])
    b = np.array([14,2,19])
    c = np.array([0,0,0,0,1,1])
    x = np.array([0,0,0,19,14,2])
    
elif test == 3: # Example for the Two-phase Method Continued (Phase 2)
    A = np.array([[3,2,0,0],
                  [2,-4,-1,0],
                  [4,3,0,1]])
    b = np.array([14,2,19])
    c = np.array([2,3,0,0])
    x = np.array([4,1,2,0])
    
elif test == 4: 
    A = np.array([[2,0,0,1,0],
                  [1,1,2,0,1]])
    b = np.array([4,2])
    c = np.array([0,0,0,1,1])
    x = np.array([0,0,0,4,2])
        
# initialize as not optimal and is not unbounded
isOptimal = False; isUnbounded = False

# continue iterating until optimal point is found
while (not isOptimal) and (not isUnbounded):
    x,isOptimal,isUnbounded = simplex_iteration(A,b,c,x)
    

# check with linprog
# specify method 
method = 'highs'

# linprog options
options = {'disp':False}

# solve the optimization problem
nx = len(c)

# bounds
LB = np.zeros((nx,))
UB = nx*[None]

bounds = np.vstack([LB,UB]); bounds = bounds.T

# solve the problem
res = linprog(c,A_eq = A,b_eq = b,bounds = bounds,method = method, options = options)

# extract result
X = res.x

print('')
print('x with linprog = {}'.format(X))
print('')
print('x with simplex = {}'.format(x))