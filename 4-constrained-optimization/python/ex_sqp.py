import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def disp_helper(name,number):
    
    # form string
    string = name + ' = {}'
    
    # display string
    print(string.format(number))
    
def plot_helper_1(F,Q):
    
    fig,ax = plt.subplots(1)
    ax.set_xlim([-1.5,0])
    ax.set_ylim([0,1.5])
    
    N = 1000
    
    x1 = np.linspace(-1.5,0,N)
    x2 = np.linspace(0,1.5,N)
    
    X1,X2 = np.meshgrid(x1,x2)
    
    # reshape
    # X1 = np.reshape(X1,[N*N,],order = 'F')
    # X2 = np.reshape(X2,[N*N,],order = 'F')
    
    F_ = F(X1,X2) 
    ax.contour(X1,X2,F_,50)
    ax.set_xlabel('x1'); ax.set_ylabel('x2')
    
    return ax
    
    
    
## create functions and serivatives
# symbolic functions and derivatives
F = lambda x1,x2: np.exp(3*x1) + np.exp(-4*x2)
Q = np.eye(2)
dF = lambda x1,x2: np.array([2*np.exp(3*x1) + -4*np.exp(-4*x2)])
d2F = lambda x1,x2: np.array([[9*np.exp(3*x1),0],[0,16*np.exp(-4*x2)]])

G = lambda x1,x2: x1**2+x2**2-1
dG = lambda x1,x2: np.array([2*x1,2*x2])

dL = lambda x1,x2,l: np.array([3*np.exp(3*x1) -2*l*x1,-4*np.exp(-4*x2)-2*l*x2])
d2L = lambda x1,x2,l: np.array([[9*np.exp(3*x1)-2*l,0],[0,16*np.exp(-4*x2)-2*l]])

## setup

test = 1 # see below

if test == 1:
    # initial point
    x = np.array([-1,1]); l = -1
elif test == 2:
    # initial point
    x = np.array([-0.75,0.1]); l = 1
elif test == 3:
     # initial point
    x = np.array([-1,1]); l = 0
    
#plot quadratic model (in 3d)
modelflag = False

# tolerences
ConstraintTolerance = 1e-10; # constraint tolerance
OptimalityTolerance = 1e-10; # optimality tolerance
MaxIterations = 100; # maximum number of iterations

## Sequential Quadratic Programming Method
# problem information
n = len(x)
m = 1

# merit penalty parameter
rho = 10

# create plot
ax=plot_helper_1(F,Q)