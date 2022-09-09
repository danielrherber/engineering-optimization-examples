import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt
from numpy.linalg import norm

def newton_iteration(xk,f,fd):

    # compute step size from linear system A*p = -b
    A = fd(xk)
    b = f(xk)

    p = solve(-A,b)

    xk1 = xk + p

    return xk1

def plot_start_1d(X,f,x_min,x_max):

    # initialize plot
    fig,ax = plt.subplots(1)
    ax.set_xlabel('x',fontsize = 16)
    ax.set_ylabel('f(x)',fontsize = 16)
    ax.set_xlim([x_min,x_max])
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    # evaluate function
    ng = 10000
    Xg = np.linspace(x_min,x_max,ng)
    Fg = np.zeros(ng)

    for i in range(ng):
        Fg[i] = f(Xg[i])

    # plot
    ax.plot(Xg,Fg,'k',linewidth = 1.5)
    ax.plot(X[0],f(X[0]),'g.',markersize = 16)
    ax.hlines(0,xmin = x_min,xmax = x_max,color = 'b')

    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    return ax

def plot_start_2d(X,f,x_min,x_max):

    # initialize plot
    fig,ax = plt.subplots(1)
    ax.set_xlabel('x_1',fontsize = 16)
    ax.set_ylabel('x_2',fontsize = 16)
    ax.set_xlim([x_min,x_max]);ax.set_ylim([x_min,x_max])

    ng = 100

    # evaluate function
    Xg = np.linspace(x_min,x_max,ng)
    [X1,X2] = np.meshgrid(Xg,Xg)

    X1 = np.reshape(X1,(ng*ng,),order = 'F')
    X2 = np.reshape(X2,(ng*ng,),order = 'F')

    Fg = np.zeros(ng*ng)


    for i in range(ng*ng):
        x_ = [X1[i],X2[i]]
        Fg[i] = norm(f(x_))

    # reshape
    Fg = np.log(Fg)
    Fg = np.reshape(Fg,(ng,ng),order = 'F')
    X1 = np.reshape(X1,(ng,ng),order = 'F')
    X2 = np.reshape(X2,(ng,ng),order = 'F')

    # plot
    ax.contourf(X1,X2,Fg)
    ax.plot(X[0],X[1],'g.',markersize = 16)

    return ax

example_number = 3

if example_number == 1: # one-dimensional example
    f = lambda x: 7*x**4 + 3*x**3 + 2*x**2  + 9*x**1 + 4
    fd = lambda x: 28*x**3 + 9*x**2 + 4*x**1 + 9
    n = 8
    X = np.zeros((n+1,)) # initial point
    X[0] = 0
    scalar_flag = True
    xmin = -1.5; xmax = 0.5 # minimum and maximum value of x for plotting

elif example_number == 2: # one dimensional example when f'(x*) = 0
    f = lambda x: x**4 - 7*x**3 + 17*x**2  -17*x + 6
    fd = lambda x: 4*x**3 -21*x**2 + 34*x -17
    n = 25
    X = np.zeros((n+1,))
    X[0] = 1.1 # initial point (try 1.1 and 2.1)
    xmin = 0.5; xmax = 3 # minimum and maximum x value for plotting
    scalar_flag = True #

elif example_number == 3:
    f = lambda x: np.array([3*x[0]*x[1] + 7*x[0] + 2*x[1] - 3,
                            5*x[0]*x[1] - 9*x[0] - 4*x[1] + 6])
    fd = lambda x: np.array([[3*x[1] + 7, 3*x[0] + 2],
                             [5*x[1]-9,5*x[0]-4]])
    scalar_flag = False
    n = 10
    X = np.zeros((n+1,2)) # initial point
    X[0,:] = [1,2]
    xmin = -6; xmax = 6 # minimum and maximum x value for plotting

elif example_number == 4:
    f  = lambda x: (np.exp(x) -np.exp(-x))/(np.exp(x) + np.exp(-x))
    fd = lambda x: 4*np.exp(x)/(np.exp(2*x) + 1)**2
    n = 10
    X = np.zeros((n+1,))
    X[0] = 1.1 # initial point
    xmin = -6; xmax = 6; # minimum and maximum x value for plotting
    scalar_flag = True

# create plot
if scalar_flag: # 1-d case
    ax = plot_start_1d(X,f,xmin,xmax)
else:
    ax = plot_start_2d(X[0,:],f,xmin,xmax)

# go through each iteration
for i in range(n):

    # compute newton iteration
    if scalar_flag:
        x_ = newton_iteration(X[i], f, fd)
        f_ = f(x_)

        # assign
        X[i+1] = x_

        # plot current point
        ax.plot(x_,f_,'r.',markersize = 16)

        # display stuff
        print('x:{} f(x):{}'.format(x_,f_))

    else:
        x_ = newton_iteration(X[i,:], f, fd)
        f_ = f(x_)

        # assign
        X[i+1,:] = x_

        # plot current point
        ax.plot(x_[0],x_[1],'r.',markersize = 16)

        # display stuff
        print('x_1:{}, x_2:{} f(x):{}'.format(x_[0],x_[1],norm(f_)))

# show plot
plt.show()