import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve,norm
from scipy.linalg import solve as solve2
# from sympy.utilities.autowrap import ufuncify
# from sympy import symbols, exp,Lambda,diff,hessian
from scipy.optimize import fminbound
from matplotlib.pyplot import quiver

def disp_helper(name,number):
    
    # form string
    string = name + ' = {}'
    
    # display string
    print(string.format(number))
    
def plotcircle(r,x,y,ax):
    
    # parameter between 0 and 2*pi
    th = np.arange(0,2*np.pi,np.pi/100)
    
    # complex-valued circle equation
    f = r*np.exp(1j*th) + x+1j*y
    
    # plot
    ax.plot(f.real,f.imag,linewidth=2)
    
    return ax

def plot_helper(flag,data):
    
    if flag == 1:
        # extract 
        F_ = data[0]
        modelflag = data[1]
        n = data[2]
        limits = data[3]
        
        # create a grid of points
        N = 600
        x1 = np.linspace(limits[0],limits[1],N)
        x2 = np.linspace(limits[2],limits[3],N)
        
        X1,X2 = np.meshgrid(x1,x2)
        
        F__ = F_(X1,X2)
        
        # create contour plot
        fig,ax = plt.subplots(1)
        ax.contourf(X1,X2,F__,50)
        
        if ~modelflag:
            ax.axis('equal')
            
        # limits
        xlim = [np.min(X1),np.max(X1)]
        ylim = [np.min(X2),np.max(X2)]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        # font size
        ax.set_xlabel('X1',fontsize = 15)
        ax.set_ylabel('X2',fontsize = 15)
        
        
    elif flag == 2:
        
        # extract
        x0 = data[0]
        k = data[1]
        D0 = data[2]
        modelflag = data[3]
        f0 = data[4]
        g0 = data[5]
        h0 = data[6]
        ax = data[7]
        
        ax.plot(x0[0],x0[1],'.',color = 'w',markersize = 20)
        
        # plot current trust-region circle
        ax = plotcircle(D0,x0[0],x0[1],ax)
        
        # only if we are visualizing the quadratic model
        if modelflag:
            
            # create grid of points around current point
            D0_ = np.linspace(-D0,D0,600)
            D0_1,D0_2 = np.meshgrid(D0_,D0_)
            
            #values of the quadratic model
            Q0 = np.zeros((600*600,))
            D0_1 = np.reshape(D0_1,[600*600,],order = 'F')
            D0_2 = np.reshape(D0_2,[600*600,],order = 'F')
            
            for idx in range(600*600):
                P0 = np.array([D0_1[idx],D0_2[idx]])
                Q0[idx] = f0 + g0.dot(P0) + 0.5*P0.dot(h0.dot(P0))
                
            D0_1 = np.reshape(D0_1,[600,600],order = 'F')
            D0_2 = np.reshape(D0_2,[600,600],order = 'F')
            Q0 = np.reshape(Q0,[600,600],order = 'F')
            
            
            ax.plot_surface(x0[0] + D0_1, x0[1] + D0_2,Q0,alpha = 0.75)
           

            
    elif flag == 3:
         
        # extract
        modelflag = data[0]
        x0 = data[1]
        q0 = data[2]
        k = data[3]
        p0 = data[4]
        pN = data[5]
        g0 = data[6]
        ax = data[7]
        
        # step direction and length
        quiver(x0[0]-p0[0],x0[1]-p0[1],p0[0],p0[1],scale = 0.9,scale_units = 'inches',linewidth=1.5,color = 'r')
        
        # newton direction
        quiver(x0[0]-p0[0],x0[1]-p0[1],pN[0]/norm(pN),pN[1]/norm(pN),scale = 1.5,color = 'b',scale_units = 'inches',linewidth = 0.2)
        
        # negative gradient direction
        quiver(x0[0]-p0[0],x0[1]-p0[1],-g0[0]/norm(g0),-g0[1]/norm(g0),color = 'g',scale = 1.5,scale_units = 'inches',linewidth = 0.2)

        # plot new point
        ax.plot(x0[0],x0[1],'.',markersize = 20)
         
         
         
        
    return ax



# example number
test = 2

if test == 1:

    f = lambda x1,x2: x1**4 + 2*x1**3 + 24*x1**2 + x2**4 + 12*x2**2
    g = lambda x1,x2: np.array([4*x1**3 + 6*x1**2 + 48*x1,4*x2**3 + 24*x2])
    h = lambda x1,x2: np.array([[12*x1**2+12*x1+48,0],[0,12*x2**2 + 24]])
    
    x0 = np.array([2,1])
    limits = [-1,3,-1,1.5]
    
    # number of iterations
    n = 4
    
elif test == 2:
    
   f = lambda x1,x2: x1**4 + 2*x1**3 + 24*x1**2 + x2**4 + 12*x2**2 + x1**6*x2**4
   g = lambda x1,x2: np.array([4*x1**3 + 6*x1**2 + 48*x1 + 6*x1**5*x2**4,4*x2**3 + 24*x2 + 4*x1**6*x2**3])
   h = lambda x1,x2: np.array([[12*x1**2+12*x1+48+30*x1**4*x2**4,24*x1**5*x2**3],
                               [24*x1**5*x2**3,12*x1**6*x2**2 + 12*x2**2 + 24]])
   
   x0 = np.array([2.5,1])
   limits = [-1,3,-1,1.5]
   
   
   # number of iterations
   n = 7 # need to set n = 7 to converge to optimal
        

# do you want to visualize the trust region model?
modelflag = False # Setting true doesnt work for now

# initial trust region size 
D0 = 1

# algorithm parameters
eta = 0.75
mu = 0.25
eps = np.finfo(float).eps

data = [f,modelflag,n,limits]

# plot stuff
ax = plot_helper(1,data)

# display initial function value
f0 = f(x0[0],x0[1])
disp_helper('f(x)',f0)

# go through each iteration
for k in range(n):
    
    print('')
    
    # current function value, gradient hessian
    f0 = f(x0[0],x0[1])
    g0 = g(x0[0],x0[1])
    h0 = h(x0[0],x0[1])
    
    data = [x0,k,D0,modelflag,f0,g0,h0,ax]
    ax = plot_helper(2,data)
    
    # compute Newton direction (solving a set of linear equations)
    pN = -solve2(h0,g0)

    # check norm of newton direction
    if norm(pN) < D0:
        
        # assign step as the newton step
        p0 = pN
        print('Using Newton Step')
        
    else:
        
        # solve for l
        obj = lambda l: (norm(solve(h0 +l*np.eye(2),-g0))-D0)**2
        
        l = fminbound(obj,0,1e10,xtol = 1e-12)
        
        # compute trust-region step
        p0 = -solve(h0 + l*np.eye(2),g0)
        
    # compute, for this step, the trust-region model value
    q0 = f0 + g0.dot(p0) + 0.5*p0.dot(h0.dot(p0))
    
    # compute function value at x0 + p0 
    f1 = f(x0[0]+p0[0],x0[1]+p0[1])
    
    # compute ratio of actual to predicted reduction
    rho0 = (f0 - f1)/(f0 - q0)
    
    # display function value and prediction ratio
    disp_helper('f(x)',f1)
    disp_helper('rho',rho0)
    
    # determine if the step was successfull or not
    if norm(f0-q0) < eps:
        print("f and q the same")
        x0 = x0 + p0
        print('Stopping')
        break
    elif rho0 <= mu:
        D0 = 0.5*D0 # decrease
        print("Trust region size decreased")
        print("Step unsucessful")
    elif rho0 >= eta:
        D0 = 2*D0 # increase
        print('Trust region size increased')
        x0 = x0+p0
        print('Step successful')
    else:
        print('Trust region size the same')
        x0 = x0+p0
        print('step successful')
        
    # plot stuff
    ax = plot_helper(3,[modelflag,x0,q0,k,p0,pN,g0,ax])


plt.show()
    
