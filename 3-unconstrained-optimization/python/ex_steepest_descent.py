import numpy as np 
import matplotlib.pyplot as plt
from numpy.linalg import norm,solve,qr
import numpy.random as random

def create_rand_pd_matrix(d,s):
    
    # set random seed
    random.seed(s)
        
    # size of the matrix 
    n = len(d)
    
    # get a random rotation matrix
    Q,r = qr(random.random((n,n)))

    # apply the transformation
    d_ = np.diag(d)
    A = Q.dot(d_.dot(Q))
    
    return A
    

test = 2

if test == 1: # example from the book
    Q = np.diag([1,5,25])
    c = np.array([-1,-1,-1])
    x0 = np.array([0,0,0])
    n = 217
    
elif test == 2: # cond(Q) is large
    d = np.array([1,25])
    s = 3525
    Q = create_rand_pd_matrix(d,s)
    c = np.array([-1,-1])
    x0 = np.array([0,0])
    n = 20
  
elif test == 3: # cond(Q) is smallish
    d = np.array([1,2])
    s = 3243
    Q = create_rand_pd_matrix(d,s)
    c = np.array([-1,-1])
    x0 = np.array([0,0])
    n = 20
    
elif test == 4: # cond(Q) = 1
    Q = np.diag([5,5])
    c = np.array([-1,-1])
    x0 = np.array([0,0])
    n = 5
    
# optimal solution 
x_opt = solve(Q,c)

# check if this is a 2d problem for plotting
if len(Q) == 2:
    n2flag = True
else:
    n2flaf = False

# setup figure
if n2flag:
    fig,ax = plt.subplots(1)
    ax.set_xlabel('x',fontsize = 16)
    ax.set_ylabel('y',fontsize = 16)
    
# assign initial point
x = x0

# initialize
f = n*[None]

# go through each iteration
for k in range(n):
    
    # store old value
    xold = x
    
    # compute current function value
    if k < 20:
        f[k] = x.dot(Q.dot(x))/2 - c.dot(x)
        
    # steepest-descent direction
    p = -(Q.dot(x) - c)
    
    # norm of the gradient
    norm_g = norm(Q.dot(x) - c)
    print('norm(g) = {}'.format(norm_g))
    
    # exact line search 
    alpha = p.dot(p)/(p.dot(Q.dot(p)))
    
    # next point 
    x = x + alpha*p 

    # plot if 2 dimensional
    if n2flag:
        
        # plot line between old and new point
        ax.plot([xold[0],x[0]],[xold[1],x[1]],'-b',linewidth = 1)
        
        # plot current point
        ax.plot(xold[0],xold[1],'.r',markersize = 12)
        
        
# plot if it is a two dimensional problem
if n2flag:
    
    # create grid
    N = 1000
    
    x = np.linspace(-1.5,0.5,N)
    y = np.linspace(-1,0.5,N)
    
    X,Y = np.meshgrid(x,y)
    
    X = X.reshape([N*N,])
    Y = Y.reshape([N*N,])
    C = np.zeros((N*N,))
    
    for k in range(N*N):
        x = np.array([X[k],Y[k]])
        C[k] = x.dot(Q.dot(x))/2 - c.dot(x)
    
    # reshape
    X = np.reshape(X,(N,N),order = 'F')
    Y = np.reshape(Y,(N,N),order = 'F')
    C = np.reshape(C,(N,N),order = 'F')
    
    ax.contour(X,Y,C)
    
    ax.plot(x_opt[0],x_opt[1],'.g',markersize = 12)
    
    plt.show()