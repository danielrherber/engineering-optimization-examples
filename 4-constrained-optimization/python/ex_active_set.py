import numpy as np
from numpy.linalg import solve,pinv,norm 
from matplotlib.patches import Polygon
from scipy.linalg import null_space as null
import matplotlib.pyplot as plt
from scipy.sparse import find

def disp_helper(name,number,n):
    
    # default value of the number of digits
    if len(n) == 0:
        n = [5]
        
    # form string
    string = name + ' = {}'
    
    # display string
    print(string.format(number))

def plot_helper(Q,c,d,A,b):
    
    fig,ax = plt.subplots(1)
    N = 100
    
    x1 = np.linspace(-0.5,4.5,N)
    x2 = np.linspace(-0.5,3,N)
    ax.set_xlim([np.min(x1),np.max(x1)])
    ax.set_ylim([np.min(x2),np.max(x2)])
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    
    # create grid
    X1,X2 = np.meshgrid(x1,x2)
    
    # reshape
    X1 = np.reshape(X1,[N*N,],order = 'F')
    X2 = np.reshape(X2,[N*N,],order = 'F')
    
    F_ = np.zeros((N*N,))
    
    for i in range(N*N):
        x_ = np.array([X1[i],X2[i]])
        
        F_[i] = 0.5*x_.dot(Q.dot(x_)) + c.dot(x_) + d
        
    # reshape
    X1 = np.reshape(X1,[N,N],order = 'F')
    X2 = np.reshape(X2,[N,N],order = 'F')
    F_ = np.reshape(F_,[N,N],order = 'F')
    
    ax.contourf(X1,X2,F_,50)
    
    # hard coded constraint boundary
    ax.plot([4/3,0,4,4/3],[8/3,0,0,8/3],color = 'g')
    
    return ax
    
    


# problem data f(x) = 1/2*x'*Q*x + c'*x + d and Ax >= b
Q = np.array([[1,0],[0,2]])
c = np.array([-3,-4])
d = 17/2 
A = np.array([[2,-1],[-1,-1],[0,1]])
b = np.array([0,-4,0])

# create the plot
ax = plot_helper(Q,c,d,A,b)

# initial feasible point 
x = np.array([0,0])
ax.plot(x[0],x[1],'.',color = 'b',markersize = 20)
ax.text(x[0]-0.2,x[1],'0',color = 'w',fontsize = 14)

# tolerences
ConstraintTolerence = 1e-14
OptimalityTolerence = 1e-14
MaxIterations = 100

# determine initial working set
Iw = find(abs(A.dot(x)-b) <= ConstraintTolerence)
Iw = Iw[1]

# update working set
Ab = A[Iw,:] # constraint matrix
Zb = np.array([0,0]) # null space of the constraint matrix 
Abr = pinv(Ab) # right inverse of the constraint matrix

# number of constraints
m = len(A)

# initial optimality flag
IsOptimal = False

for iter_ in range(MaxIterations):
    # --- optimality test
    # compute gradient 
    g = Q.dot(x) + c
    
    matrixflag = False
    # check reduced gradient
    while norm(np.dot(Zb,g)) <= OptimalityTolerence :
        
        # compute (approximate) Lagrange multipliers
        Lb = np.dot(Abr.T,g)
        
        # find all negative multipliers
        In = Lb < 0

        # stop if all multipliers are nonnegative
        if not any(In):
            print('Optimal')
            IsOptimal = True 
            break 
        
        # select most negative multiplier 
        Ir = np.argmin(Lb[In])
        
        # remove the constraint from the working set
        Iw = np.delete(Iw,Ir)
        
        # update working set
        Ab = A[Iw,:] # constraint matrix
        
        
        if len(Ab) == 0:
            Zb = np.eye(Ab.shape[1])
            matrixflag = True      
        else:
            Zb = np.squeeze(null(Ab)) # null space of the constraint matrix 
             
            
        Abr = pinv(Ab) # right inverse of the constraint matrix

    # check if optimal
    if IsOptimal:
        break

    #---The search Direction
    # compute reduced Newton search direction
    if matrixflag:
        p_ = solve(Zb.dot(Q.dot(Zb)),Zb.dot(g))
        p = -Zb.dot(p_)
        
    else:
        p_ = Zb.dot(g)/Zb.dot(Q.dot(Zb))
        p = -Zb*p_
   
    # --- The STep
    # initialize maximum step as the the newton step
    alpha = 1
    
    # go through each constraint
    for k in range(m):
        
        # only if the constraint is not in the working set
        if not any(k == Iw):
            
            # extract current constraint row
            a = A[k,:]

            # check if this is a descent direction
            if a.dot(p) < 0:
                
                # compute maximum step length with ratio test
                alpha_ = -(a.dot(x) - b[k])/(a.dot(p))
                
                # update maximum step length
                alpha = np.min([alpha_,alpha])
                
    #---The Update
    # take step 
    x = x + alpha*p
    
    # update working set (this add all approximately active constraints)
    Iw = find(abs(A.dot(x)-b) <= ConstraintTolerence)
    Iw = Iw[1]
    
    # update working set
    Ab = A[Iw,:] # constraint matrix
    Zb = np.squeeze(null(Ab)) # null space of the constraint matrix 
    Abr = pinv(Ab) # right inverse of the constraint matrix
    
    # display stuff
    disp_helper("---Iteration",iter_,[])
    disp_helper("x",x,[])
    disp_helper("Iw",Iw,[])
    disp_helper("Zb",Zb,[])
    ax.plot(x[0],x[1],'.r',markersize=24)
    ax.text(x[0]-0.2,x[1],str(iter_+1),color = 'w',fontsize = 14)
    
    
    # compute Lagrange multipliers
    g = Q.dot(x) + c
    Lb = np.dot(Abr.T,g)
    E = Ab*Lb - g
    disp_helper("Lb",Lb,[])
    disp_helper("Lb Error",E,[])
    print('')
    
    if norm(E) > 1e-8:
        ax.text(x[0]+0.15,x[1],'inconsistent multipliers',color = 'w',fontsize = 14)
    
    plt.show()