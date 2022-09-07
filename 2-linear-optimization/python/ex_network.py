import numpy as np
from scipy.optimize import linprog
import networkx as nx # needed for visualization
import matplotlib.pyplot as plt

# problem data
c = np.array([12,15,9,8,7,6,8,4,5,3,11])

# initialize
Aeq = np.zeros((7,11))
beq = np.array([50,0,0,0,0,-20,-30])

Aeq[0,0:2] = [1,1]
Aeq[1,0:4] = [-1,0,1,1]
Aeq[2,1:6] = [-1,0,0,1,1]
Aeq[3,2:8] = [-1,0,-1,0,1,1]
Aeq[4,5:10] = [-1,-1,0,1,1]
Aeq[5,3:] = [-1,0,0,0,-1,-1,0,1]
Aeq[6,9:] = [-1,-1]

LB = np.zeros((11,))
UB = 30*np.ones((11,))

# combine bounds and reshape
bounds = np.vstack([LB,UB]); bounds = bounds.T

# specify method
method = 'highs'

# linprog options
options = {'disp':True}

# solve the optimization problem
res = linprog(c,A_eq = Aeq,b_eq = beq,bounds = bounds,method = method, options = options)


# create directed graph
G = nx.DiGraph()
G.add_node('1')
G.add_node('2')
G.add_node('3')
G.add_node('4')
G.add_node('5')
G.add_node('6')
G.add_node('7')

# edge list
edge_list = [('1','2'), ('1','3'), ('2','4'),('2','6'),('3','4'),('3','5'),
             ('4','5'),('4','6'),('5','6'),('5','7'),('6','7') ]

# edge label
edge_labels = dict(zip(edge_list,list(c)))

G.add_edges_from(edge_list)

pos = nx.kamada_kawai_layout(G)
nx.draw(G,pos = pos,with_labels = True)
nx.draw_networkx_edges(
    G,
    pos,
    edgelist = [('1','2'),('1','3'),('2','6'),('3','5'),('5','7')],
    width = 4,
    edge_color = 'tab:red')

nx.draw_networkx_edge_labels(
    G,
    pos,
    edge_labels = edge_labels
    )

plt.show()