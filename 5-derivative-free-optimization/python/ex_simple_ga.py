import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.animation import FuncAnimation 
from numpy.random import seed
from scipy.optimize import minimize,NonlinearConstraint,fsolve
from scipy.sparse import find

def roulette_wheel_selection(f):
    
    # determine number of points in a generation
    np_ =  len(f)
    np__ = int(np_/2)
    
    # minimimum, maximum , and spread of the function values
    fmin = np.min(f)
    fmax = np.max(f)
    Fdelta = 1.1*fmax -0.1*fmin 
    
    # scaled fitness values
    F = (-f+Fdelta)/(np.max([1,Fdelta-fmin]))
    
    # normalized cumulative sum of the scaled fitness values
    S = np.cumsum(F)/np.sum(F)

    # random matrix for creating the pairs
    R = np.random.random([np__,2])
    
    # determine the index of the parent
    parent_indices = np.zeros([np__,2])
    
    for i in range(np__):
        
        # determine the first parent
        ind1 = find(S >= R[i,0]); ind1 = ind1[1]; ind1 = ind1[0]
        
        # determine the second parent
        ind2 = find(S >= R[i,1]); ind2 = ind2[1]; ind2 = ind2[0]

        parent_indices[i,0] = ind1
        parent_indices[i,1] = ind2
    
    return parent_indices.astype(int)

def linear_crossover(f,x,parent_indices,alpha):
    
    # determine number of the points in a generation
    np_ = len(f)
    np__ = int(np_/2)
    
    # determine children
    children = np.zeros((np_,2))
    
    # go through each pair of parents
    for i in range(np__):
        
        # extract parent indices
        parent1_index = parent_indices[i,0]
        parent2_index = parent_indices[i,1]

        # flip points if ordering is wrong
        if f[parent2_index] > f[parent1_index]:
            temp = parent1_index
            parent1_index = parent2_index
            parent2_index = temp
            
        # extract parents
        xp1 = x[parent1_index,:]
        xp2 = x[parent2_index,:]
        
        # child 1
        children[2*i-1,:] = 0.5*xp1+0.5*xp2
        
        # child 2
        children[2*i,:] = (1+alpha)*xp2 -alpha*xp1
        
    return children

# implementation of uniform random mutation on 2 variables
def uniform_random_mutation(children,p,Delta):

    # determine number of points in a generation
    np_ = len(children)

    #  random [0,1] matrix for determining if mutation should occur
    pr = np.random.random([np_,2])
    
    # random [0,1] matrix for determining how much mutation should occur
    rr = np.random.random([np_,2])
    
    # mutate some children
    children = children + (pr <=p)*(rr-0.5)*Delta 
    
    return children

# function to make it easier to display things in the command window
def disp_helper(name,number):
    
    # form string
    string = name + ' = {}'
    
    # display string
    print(string.format(number))
    
def disp_iteration(f,x,X_best,F_best,F_mean,iter_):
    
    # get best point
    f_best = np.min(f); I_best = np.argmin(f)
    
    x_best = x[I_best,:]
    f_mean = np.mean(f)
    
    X_best[iter_,:] = x_best 
    F_best[iter_] = f_best 
    F_mean[iter_] = f_mean
    
    # display stuff
    disp_helper("*iteration",iter_)
    disp_helper("mean(f)",f_mean)
    disp_helper("min(f)",f_best)
    
    return X_best,F_best,F_mean

def plot_helper1(obj,pen,con,lb,ub,x,max_iterations):
    
    # create grid of evalutation points
    N = 200
    x1 = np.linspace(lb[0],ub[0],N)
    x2 = np.linspace(lb[1],ub[1],N)
    X1,X2 = np.meshgrid(x1,x2)
    
    # evaluate the function
    F_ = obj(X1,X2)
    
    if len(pen) != 0:
        pen = pen[0]
        F_ = F_ - pen(X1,X2)
        
    fig,ax = plt.subplots(1)
    ax.contour(X1,X2,F_,30)
    ax.set_xlim([lb[0],ub[0]]); ax.set_ylim([lb[1],ub[1]])
    ax.set_xlabel('x1',fontsize=18); ax.set_ylabel('x2',fontsize = 18)
    
    # plot initial population
    ax.plot(x[:,0],x[:,1],'g.',markersize = 10)
    
    # plot colors
    parent_color = [0.5,0.5,0.5]
    
    
    return fig,ax,parent_color

#-------Example starts from here----------------------    

test = 3

if test == 1:
    fun = 1
    max_iterations = 25 # maximum number of iterations
    np_ = 30 # population size( even number for simplicity)
    seed(876044) # set random seed for repeatable results
    
elif test == 2:
    fun = 2
    max_iterations = 30 # maximum number of iterations
    np_ = 30 # population size (even number for simplicity)
    seed(876044) # set random seed for repeatable results
    
elif test == 3:
    fun = 3
    max_iterations = 30 # maximum number of iterations
    np_ = 30 # population size (even number for simplicity)
    seed(876044) # set random seed for repeatable results
    
    
# algorithm parameters
p = 0.01 # mutation probability
Delta = 0.5 # maximum mutation amount
alpha = 0.5 # linear crossover parameter

#create objective function and derivatives
pen = []; con = [];

# define objective function
if fun == 1:
    
    # spring examples
    l1 = 12; l2 = 8; k1 = 1; k2 = 10; mg = 7; # problem data
    obj = lambda x,y: k1*0.5*(np.sqrt((l1 + x)**2 + y**2)-l1)**2 + k2*0.5*(np.sqrt((l2-x)**2 +y**2)-l2)**2 -mg*y
    ub = np.array([16,12])
    lb = np.array([-6,-9])
   
elif fun == 2:
        
    # Jones function
    fs = lambda x,y: x**4 + y**4 -4*x**3 -3*y**3 +2*x**2 + 2*x*y; obj = fs
    ub = np.array([4,3])
    lb = np.array([-2,-2])
    
elif fun == 3:
    
    # spring example with linear constraint
    # spring examples
    l1 = 12; l2 = 8; k1 = 1; k2 = 10; mg = 7; # problem data   
    gs = lambda x,y: 2*x + y-8;
    pen = lambda x,y: 100*(2*x+y-8)**2; pen = [pen]
    fs = lambda x,y: k1*0.5*(np.sqrt((l1 + x)**2 + y**2)-l1)**2 + k2*0.5*(np.sqrt((l2-x)**2 +y**2)-l2)**2 -mg*y + 100*(2*x+y-8)**2
    ub = np.array([16,12]); obj = fs
    lb = np.array([-6,-9])
    
    
    
# ----random initial population
x1_ = lb[0] + np.random.random((np_,))*(ub[0]-lb[0])
x2_ = lb[1] + np.random.random((np_,))*(ub[1]-lb[1])
x = np.vstack([x1_,x2_]); x = x.T

# create plot
fig,ax,parent_color = plot_helper1(obj, pen, con, lb, ub, x, max_iterations)


# get optimal solution using scipy minimize so that it can be plotted
if fun == 1:
    obj_ = lambda x: obj(x[0],x[1])
    res = minimize(obj_,x0 = [1,1],method='Nelder-Mead',options={'disp': False})
    X_optimal = res.x
    
    # plot optimal point
    ax.plot(X_optimal[0],X_optimal[1],'.',markersize = 30,color = 'tab:pink')
    
elif fun == 2:
    obj_ = lambda x: obj(x[0],x[1])
    res = minimize(obj_,x0 = [2.5,-1],method='Nelder-Mead')
    X_optimal = np.array([[-0.4495,2.2928],[2.4239,1.9219]])
    X_optimal = np.vstack([res.x,X_optimal])
    
    # plot optimal point
    ax.plot(X_optimal[:,0],X_optimal[:,1],'.',markersize = 30,color = 'tab:pink')
    
elif fun == 3:   
    obj_ = lambda x: obj(x[0],x[1])
    con = lambda x:gs(x[0],x[1]); nlc = NonlinearConstraint(con,0,0)
    res = minimize(obj_, [0,0], method='SLSQP', constraints=nlc)
    
    # extract
    X_optimal = res.x
    
    # plot optimal point
    ax.plot(X_optimal[0],X_optimal[1],'.',markersize = 30,color = 'tab:pink')
    
    # solve for constraints
    con1 = lambda x: gs(lb[0],x)
    y1 = fsolve(con1,lb[1])
    
    con2 = lambda x: gs(ub[0],x)
    y2 = fsolve(con2,lb[1])
    
    ax.plot([lb[0],ub[0]],[y1,y2],'-',linewidth = 2, color = 'tab:green')
    
       
# initialize
X_best = 0*np.ones((max_iterations,2)); F_best = 0*np.ones((max_iterations,)); F_mean = 0*np.ones((max_iterations,))   
  
# go through each generation
def ga_iteration(iter_):
    
    global x,X_best,F_best,F_mean,children
    
    # evaluate population
    f = obj(x[:,0],x[:,1])
    
    # display and calulate things related to this iteration
    X_best,F_best,F_mean = disp_iteration(f,x,X_best,F_best,F_mean,iter_);
    
    # selection step
    parent_indices = roulette_wheel_selection(f)
    
    # crossover step
    children = linear_crossover(f,x,parent_indices,alpha)
    
    # mutation step
    children = uniform_random_mutation(children,p,Delta)
    
    # fix points outside the bounds
    children[(children[:,0]<lb[0]),0] = lb[0]
    children[(children[:,1]<lb[1]),1] = lb[1]
    
    children[(children[:,0]>ub[0]),0] = ub[0]
    children[(children[:,1]>ub[1]),1] = ub[1]
    
    # set current children are the next population
    x = children
    
    
def animate(i):

    title = 'Iteration {:02d}'.format(i)
    
    ax.plot(x[:,0],x[:,1],'.',markersize = 10,color = 'tab:grey')
    
    # golden section update
    ga_iteration(i)

    # set picture
    ax.set_title(title)
    
    # plot interior points
    ax.plot(children[:,0],children[:,1],'.',markersize = 10,color = 'tab:red')
    
    
    return ax


# create animation and save
savename = 'genetic_algorithm' + '_ex_' + str(fun) + '.gif'
anim = FuncAnimation(fig,animate,frames = list(range(max_iterations)),interval=500,blit=False,repeat=False)
anim.save(savename, dpi=120, writer="imagemagick")

# plot convergence behavior
fig2,ax2 = plt.subplots(1)
ax2.plot(np.arange(0,max_iterations),F_mean,'b.-',linewidth = 2,label = 'Best')
ax2.plot(np.arange(0,max_iterations),F_best,'r.-',linewidth = 2,label = 'Mean')
ax2.legend()
ax2.set_xlabel('Iteration Number',fontsize=18); ax2.set_ylabel('f',fontsize = 18)
fig2.show()