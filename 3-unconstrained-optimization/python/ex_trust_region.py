import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve
from sympy.utilities.autowrap import ufuncify
from sympy import symbols, exp,Lambda,diff,hessian

test = 1

if test == 1:

    x1,x2 = symbols('x1 x2')
    f = x1**4 + 2*x1**3 + 24*x1**2 + x2**4 + 12*x2**2
    g = diff(f,[x1,x2])
    
    x0 = np.array([2,1])
    limits = [-1,3,-1,1.5]
    
    
