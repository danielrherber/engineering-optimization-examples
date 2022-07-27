import numpy as np 
import numpy.random as random
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# get activation function parameters from optimization variables
def get_parameters(x,N):
    
    w = x[0:N] # weights
    a = x[N:2*N] # amplitude
    b = x[2*N:] # bias
    
    
    return w,a,b

def activation(x_data,w,a,b):
    
    # reshape arrays
    w = w.reshape([-1,1]); a = a.reshape([-1,1]); b = b.reshape([-1,1])
    x_data = x_data.reshape([1,len(x_data)])
    
    # prediction
    y_pred = np.sum(w/(1 + np.exp(np.dot(a,x_data) + b)),axis = 0)
    
    return y_pred

# objective function
def objective(x,x_data,y_data,N):
    
    w,a,b = get_parameters(x,N)
    
    # predict
    y_pred = activation(x_data,w,a,b)
    
    # calculate error
    f = np.sum((y_pred-y_data)**2)/len(y_pred)
    
    return f


# set random seet
seed = 457568
random.seed(seed = seed)

n = 100
x_data = np.linspace(-10,10,n)
y_data = 0.1*x_data*np.cos(x_data) + 0.1*random.random(n)

# number of activation functions
N = 6

# starting point
X0 = random.random((3*N,))

# options
options = {'maxiter':100000000000}

# solve the unconstrained optimization problem
res = minimize(objective, X0, method='BFGS',args = (x_data,y_data,N),
               options=options)

# extract results
x = res.x

# extract optimal parameter values
w,a,b = get_parameters(x,N)

# get predictions
y_pred = activation(x_data,w,a,b)

# plot
fig,ax = plt.subplots(1)
ax.plot(x_data,y_data,'.')
ax.plot(x_data,y_pred)
plt.show()