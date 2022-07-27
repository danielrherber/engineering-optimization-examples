import numpy as np 
from scipy.optimize import linprog
import networkx as nx
import matplotlib.pyplot as plt

# problem data 
c = np.array([11,5,2,15,12,8,3,1,10])

# equality constraints
Aeq = np.array([[1,1,1,0,0,0,0,0,0],
    [0,0,0,1,1,1,0,0,0],
    [0,0,0,0,0,0,1,1,1],
    [-1,0,0,-1,0,0,-1,0,0],
    [0,-1,0,0,-1,0,0,-1,0],
    [0,0,-1,0,0,-1,0,0,-1]])
beq = np.array([1,1,1,-1,-1,-1])

# bounds
LB = np.zeros((len(c),))
UB = np.ones((len(c),))

# combine bounds and reshape
bounds = np.vstack([LB,UB]); bounds = bounds.T

# specify method 
method = 'highs'

# linprog options
options = {'disp':True}

# solve the optimization problem
res = linprog(-c,A_eq = Aeq,b_eq = beq,bounds = bounds,method = method, options = options)

# create directed graph
G = nx.DiGraph()
G.add_node('Person 1')
G.add_node('Person 2')
G.add_node('Person 3')

G.add_node('Personal Manager')
G.add_node('Budget Director')
G.add_node('Accountant')

node_list = ['Person 1','Person 2','Person 3']

edge_list = [('Person 1','Personal Manager'), ('Person 1','Budget Director'), ('Person 1','Accountant'),
                  ('Person 2','Personal Manager'), ('Person 2','Budget Director'), ('Person 2','Accountant'),
                  ('Person 3','Personal Manager'),('Person 3','Budget Director'),('Person 3','Accountant')]

G.add_edges_from(edge_list)

# create a circular layout
pos = nx.circular_layout(G)
edge_labels = dict(zip(edge_list,list(c)))

# plot the graph
nx.draw(G,pos = pos,with_labels = True)

nx.draw_networkx_edges(
    G,
    pos,
    edgelist = [('Person 3','Personal Manager'),('Person 2','Budget Director'),('Person 1','Accountant')],
    width = 3,
    alpha = 0.5,
    edge_color = 'tab:red')

nx.draw_networkx_edge_labels(
    G,
    pos,
    edge_labels = edge_labels
    )

plt.show()
