import numpy as np
import numpy.random as random
from numpy.linalg import norm,solve


# function to make it easier to display things in the command window
def disp_helper(name,number):
    
    # form string
    string = name + ' = {}'  
    
    # display string
    print(string.format(number))



test = 3
updateformulaflag = 1

if test == 1:
    Q = np.diag([2,3,4]) # x'*Q*x/2
    c = np.array([-8,-9,-8]) # c'x
    B = np.eye(len(Q)) # initial approximate hessian matrix
    x = np.zeros((len(Q),)) # initial point
    
elif test == 2:
    random.seed(3253) # set random seed
    q = random.random((3,3))
    Q = q + q.T # x'*Q*x/2
    c = np.array([-8,-9,-8]) # c'x
    B = np.eye(len(Q)) # initial approximate hessian matrix
    x = np.zeros((len(Q),)) # initial point
    
elif test == 3:
    random.seed(786796)
    n = 20
    q = random.random((n,n))
    Q = q + q.T + 8*np.eye(n) # x'*Q*x/2
    c = random.random(n) # c'x
    B = np.eye(n) # initial approximate hessian matrix
    x = np.zeros((n,)) # initial point
    
    
    
# exact gradient
G = lambda x: Q.dot(x) - c

# exact hessian
H = lambda x: Q

# maximum number of iterations
max_iterations = 100

# compute initial gradient
g = G(x)

# calculate eps
eps = np.finfo(float).eps

# go through each iteration
for i in range(max_iterations):
    
    # (i) check if xk is optional
    if norm(g,ord = np.inf) <= 100*eps:
        print('Optimal!')
        break
    else:
        disp_helper('iteration',i+1)
        disp_helper('norm(g)',norm(g,ord = np.inf))
        
    # (ii) solve for search direction
    p = -solve(B,g)
    disp_helper('p',p)

    # (iii) using exact line search 
    alpha = -(p.dot(g))/(p.dot(Q.dot(p)))
    disp_helper('alpha',alpha)
    
    x_new = x + alpha*p 
    disp_helper('x',x_new)
    
    # (iv) compute intermediate quantities
    g_new = G(x_new) # gradient
    s = x_new - x # step 
    y = g_new - g # gradient difference
    
    disp_helper('g',g_new)
    disp_helper('s',s)
    disp_helper('y',y)
    
    # (v) use symmetric rank-one update formula 
    if updateformulaflag == 1: #use symmetric rank-one update formula
        z = y - B.dot(s)
        z_ = z.reshape([-1,1])
        B_new = B + z_.dot(z_.T)/(z.dot(s))

    elif updateformulaflag == 2: # BFGS update formula
        z = B.dot(s)
        z_ = z.reshape([-1,1])
        y_ = y.reshape([-1,1])
        B_new = B - (z_.dot(z_.T))/(z.dot(s)) + (y_.dot(y_.T))/(y.dot(s))
        
    disp_helper('B',B_new)
    
    # update stored variables
    B = B_new
    g = g_new 
    x = x_new
    
    print('')
    

# compute error compared to exact solution
e_norm = norm(x - solve(Q,c),ord = np.inf)
disp_helper('max abs error',e_norm)

# compute error comprared to exact hessian
e_Q = norm(B-Q,ord = np.inf)
disp_helper('error Q',e_Q)            
    
    