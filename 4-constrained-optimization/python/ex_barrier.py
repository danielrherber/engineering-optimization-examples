import numpy as np 
from scipy.optimize import minimize

def disp_helper(name,number,n):
    
    # default value of the number of digits
    if len(n) == 0:
        n = [5]
        
    # form string
    string = name + ' = {}'
    
    # display string
    print(string.format(number))
    

# problem functions
f = lambda x: x[0] -2*x[1]
g = lambda x: np.array([1 + x[0]-x[1]**2,x[1]])

# initial strictly feasible point
x = np.array([0.5,0.5])

# initial barrier parameter
mu = 1 

# barrier term
phi = lambda x: -np.sum(np.log(g(x)))

# barrier function
beta = lambda x,mu: f(x) + mu*phi(x)

# problem options
options={'xatol':0, 'fatol':0,'disp': False}

# go through each iteration
for i in range(9):
    
    res = minimize(beta,args = (mu),x0 = x,method = 'Nelder-Mead',options = options)
    
    # extract results
    x = res.x
    
    disp_helper('---iteration',i,[])
    disp_helper('mu',mu,[10])
    disp_helper('x',x,[10])
    print('')
    
    # update mu
    mu = mu/10

