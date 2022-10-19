import numpy as np
from truss import truss
from scipy.optimize import (BFGS, SR1, Bounds, NonlinearConstraint, minimize)
import matplotlib.pyplot as plt
import networkx as nx

def disp_helper(name,number,n):
    
    # default value of the number of digits
    if len(n) == 0:
        n = [5]
        
    # form string
    string = name + ' = {}'
    
    # display string
    print(string.format(number))
    
    
class problem_wrapper():
    
    '''
    For this problem it is necessary to pass extra arguments ot eh constraint function
    Passing arguments to scipy minimize for constraints is not possible. 
    Thus its necessary to wrap the constraints and objeective using a wrapper class.
    
    Note: This method works with other external solvers like pyoptsparse too
    
    '''
    
    def __init__(self,mass0,stress_limit):
        
        # assign
        self.mass0 = mass0
        self.stress_limit = stress_limit
        
    # objective function
    def objective(self,x):
        
        # comute mass
        mass,st = truss(x)
        
        # assign objective value
        f = mass/self.mass0 # scaled

        return f
    
    def constraints(self,x):
        
        # extract
        stress_limit = self.stress_limit
        
        # compute stress
        m,stress = truss(x)
        
        # inequality constraint
        c = np.hstack([stress/stress_limit-1,-stress/stress_limit -1])

        return c

# problem data
stress_limit = 25000
min_area = 0.1

# compute minimum mass for scaling
mass0,stress0 = truss(min_area*np.ones((10,)))

# call wrapper class
p = problem_wrapper(mass0,stress_limit)

# problem function and constraints
X0 = np.ones((10,))
nonlinear_constraint = NonlinearConstraint(p.constraints,-np.inf,0,jac='2-point', hess=BFGS())
LB = min_area*np.ones((10,))
UB = np.array(10*[np.inf])
bounds = Bounds(LB,UB)

res = minimize(p.objective,
                x0=X0,
                method='trust-constr',
                jac="2-point",
                hess=SR1(),
                constraints=[nonlinear_constraint],
                options={'verbose': 1},
                bounds=bounds)

# extract result
X = res.x
X_ = np.round(X,decimals=1)

# calculate the objective and constraints at minimizer
mass,stress = truss(X)

# display solution
print('')
disp_helper('X',X,[4])
disp_helper('stress',stress,[4])


# create directed graph
G = nx.DiGraph()
G.add_node('1')
G.add_node('2')
G.add_node('3')
G.add_node('4')
G.add_node('5')
G.add_node('6')


# edge list
edge_list = [('5','3'), ('3','1'), ('6','4'),('4','2'),('4','3'),('2','1'),
             ('5','4'),('6','3'),('3','2'),('4','1')]


G.add_edges_from(edge_list)

edge_labels = dict(zip(edge_list,list(X_)))

G.add_edges_from(edge_list)

pos = nx.bipartite_layout(G,nodes = ['5','3','1'],align = 'horizontal')
nx.draw(G,pos = pos,with_labels = True)


nx.draw_networkx_edge_labels(
    G,
    pos,
    edge_labels = edge_labels
    )

plt.show()