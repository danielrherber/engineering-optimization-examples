import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.animation import FuncAnimation 
from numpy.random import seed
from scipy.optimize import minimize,NonlinearConstraint,fsolve
from scipy.sparse import find
from matplotlib.pyplot import quiver


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

def plot_helper1(obj,pen,con,lb,ub,x,p,max_iterations):
    
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
    Q = quiver(x[:,0],x[:,1],p[:,0],p[:,1],color = 'g')
    
    # plot colors
    parent_color = [0.5,0.5,0.5]
    
    
    return fig,ax,Q,parent_color

def step_norm_limiter(p,p_norm_max):
    
    # compute norm of each step
    p_norm = np.sqrt(np.sum(p**2,1))
    
    # get indices that are larger than the maximum
    I = p_norm > p_norm_max

    # normalize
    p_limited = p.copy()
    p_limited[:,0] = p_norm_max*p_limited[:,0]/p_norm
    p_limited[:,1] = p_norm_max*p_limited[:,1]/p_norm

    # modify only the large steps
    p[I,:] = p_limited[I,:]
    
    return p

#-------Example starts from here----------------------    

test = 1

if test == 1:
    fun = 1
    max_iterations = 40 # maximum number of iterations
    np_ = 20 # population size( even number for simplicity)
    p_norm_max = 2 # maximum step size norm
    seed(2153) # set random seed for repeatable results
    
elif test == 2:
    fun = 2
    max_iterations = 40 # maximum number of iterations
    np_ = 10 # population size (even number for simplicity)
    p_norm_max = 0.5 # maximum step size norm
    seed(1561) # set random seed for repeatable results
    
elif test == 3:
    fun = 3
    max_iterations = 30 # maximum number of iterations
    np_ = 30 # population size (even number for simplicity)
    p_norm_max = 2 # maximum step size norm
    seed(46346) # set random seed for repeatable results
    
    
# algorithm parameters
alpha = 0.7 # inertia parameter
beta_max = 1 # self influence parameter
gamma_max = 1 # social influence parameter

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
nsamples = round(np.sqrt(np_))
x1_ = np.linspace(lb[0],ub[0],nsamples)
x2_ = np.linspace(lb[1],ub[1],nsamples)

x1_,x2_ = np.meshgrid(x1_,x2_)

x1_ = np.reshape(x1_,[nsamples**2,],order = 'F')
x2_ = np.reshape(x2_,[nsamples**2,],order = 'F')

x = np.vstack([x1_,x2_]); x = x.T;

# update population size
np_ = len(x)

# randomize initial directions (velocities)
p = 0.1*(1-2*np.random.random((np_,2)))*(ub-lb)

# limit step based on maximum step
p = step_norm_limiter(p,p_norm_max)

# create plot
fig,ax,Q,parent_color = plot_helper1(obj, pen, con, lb, ub, x, p, max_iterations)
n = 2
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
x_best_particles = np.inf*np.ones((np_,n)); f_best_particles = np.inf*np.ones(np_)
f_best_overall = np.inf; x_best_overall = [] 

def pso_iteration(iter_):
    
    global x,p,X_best,F_best,F_mean,x_best_particles,f_best_particles,f_best_overall,x_best_overall,Q
    
    # evaluate objective function
    f = obj(x[:,0],x[:,1])
    
    # display and calulate things related to this iteration
    X_best,F_best,F_mean = disp_iteration(f,x,X_best,F_best,F_mean,iter_);
    
    # update best individual points
    for i in range(np_):
        if f_best_particles[i] > f[i]:
            f_best_particles[i] = f[i]
            x_best_particles[i,:] = x[i,:]
            
    # update best swarm point
    F_best_ = np.min(f_best_particles)
    I = np.argmin(f_best_particles)
    
    if f_best_overall > F_best_:
        f_best_overall = F_best_
        x_best_overall = x_best_particles[I,:]
        
    # compute each particles step
    beta = beta_max*np.random.random((np_,))
    gamma = gamma_max*np.random.random((np_,))

    x_best_diff = (x_best_particles - x); x_overall_diff = (x_best_overall - x)
     
    p[:,0] = alpha*p[:,0]+ beta*x_best_diff[:,0] + gamma*x_overall_diff[:,0]
    p[:,1] = alpha*p[:,1] + beta*x_best_diff[:,1] + gamma*x_overall_diff[:,1]

    # limit step based on maximum step
    p = step_norm_limiter(p,p_norm_max)

    # update each particles position
    x = x + p
    
    # plot interior points
    ax.plot(x[:,0],x[:,1],'.',markersize = 10,color = 'tab:red')
    Q.remove()

    Q = ax.quiver(x[:,0],x[:,1],p[:,0],p[:,1],color = 'tab:red')
    
def animate(i):

    title = 'Iteration {:02d}'.format(i)
    
    ax.plot(x[:,0],x[:,1],'.',markersize = 10,color = 'tab:grey')
    
    # golden section update
    pso_iteration(i)

    # set picture
    ax.set_title(title)
    
    
    
    return ax
        
# create animation and save
savename = 'pso_algorithm' + '_ex_' + str(fun) + '.gif'
anim = FuncAnimation(fig,animate,frames = list(range(max_iterations)),interval=500,blit=False,repeat=False)
anim.save(savename, dpi=120, writer="imagemagick")

# plot convergence behavior
fig2,ax2 = plt.subplots(1)
ax2.plot(np.arange(0,max_iterations),F_mean,'b.-',linewidth = 2,label = 'Best')
ax2.plot(np.arange(0,max_iterations),F_best,'r.-',linewidth = 2,label = 'Mean')
ax2.legend()
ax2.set_xlabel('Iteration Number',fontsize=18); ax2.set_ylabel('f',fontsize = 18)
fig2.show() 