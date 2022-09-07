import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve
from sympy.utilities.autowrap import ufuncify
from sympy import symbols, exp,Lambda,diff,hessian

def newton_step(x,G,H):

    H_ = H(x)
    G_ = G(x)

    # solve Newton equations
    p = -G_/H_

    return p

def eval_q(x,P,q):

    # initialize
    Q = np.zeros((len(P),))

    # go through each point
    for k in range(len(P)):
        Q[k] = q(x,P[k]) # evaluate the quadratic function

    return Q

def quad_approx(F,G,H,x,p):

    F_ = float(F(x))
    G_ = float(G(x))
    H_ = float(H(x)[0])

    q = F_ + G_*p + 0.5*p*H_*p
    return q

test = 3

if test == 1:
    xk = 1 # starting point
    xlimits = np.array([-2.5,2.5]); ylimits = np.array([-0.1,3]) # plotting points
    niter = 4 # number of iterations

elif test == 2:
    xk = -1.9 # starting point
    xlimits = np.array([-3,3]); ylimits = np.array([-0.3,6]) # plotting points
    niter = 5 # number of iterations

elif test == 3:
    xk = -3 # starting point
    xlimits = np.array([-6,-0.5]); ylimits = np.array([-0.05,0.5]) # plotting points
    niter = 3 # number of iterations

# symbolic stuff (function,gradient,hessian)
x = symbols('x')
f = exp(0.5*x-1)*(x+1)**2
g = diff(f)
h = diff(g)

# create anonymous functions
F = ufuncify(x,f)
G = ufuncify(x,g)
H = ufuncify(x,h)

q = lambda x,p: F(x) + G(x)*p + 1/2*p*H(x)*p

# plot setup
N = 1000
xmin = -10; xmax = 10
X = np.linspace(xmin,xmax,N)
pmin = -1; pmax = 1 # plimits

# create plot
fig,ax = plt.subplots(1)

ax.set_xlabel('x',fontsize = 15)
ax.set_ylabel('f(x)',fontsize = 15)
ax.set_xlim(xlimits); ax.set_ylim(ylimits)

ax.plot(X,F(X),'k-',linewidth=3)
ax.plot(xk,F(xk),'r.',markersize = 30)

for i in range(niter):

    # compute step length
    p = newton_step(xk,G,H)

    # plotting
    P = np.linspace(np.min([0,p])+pmin,np.max([0,p])+pmax,N)

    ax.plot(xk + P,eval_q(xk,P,q),'-',linewidth = 2)

    # take step
    xk = xk+p

    # plot
    ax.plot(xk,q(xk-p,p),'.',markersize = 20)
    ax.plot([xk,xk],[q(xk-p,p),F(xk)],'-',linewidth=1)
    ax.plot(xk,F(xk),'.',markersize=20)

# display final point
print('x = {}'.format(xk))
plt.show()