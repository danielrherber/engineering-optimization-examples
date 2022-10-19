import numpy as np
from numpy.linalg import solve

def bar(E,A,L,phi):
    
    '''
    Compute the stiffness and stress matrix for one element

    Parameters
    ----------
    E : float
        modulus of elasticity
    A : float
        cross-sectional area
    L : float
        length of element
    phi : float
        orientation of element
    
    Outputs
    -------
    K : 4 x 4 ndarray
        stiffness matrix
    S : 1 x 4 ndarray
        stress matrix
    
    '''
    # rename
    c = np.cos(phi)
    s = np.sin(phi)
    
    # stiffness matrix
    k0 = np.array([[c**2,c*s],[c*s,s**2]])
    k0_ = np.hstack([k0,-k0])
    
    K = E*A/L*np.vstack([k0_,-k0_])
    
    # stress matrix
    S = E/L*np.array([-c,-s,c,s])

    return K,S

def node2idx(node,DOF):
    
    '''
    Computes the appropriate indices in the global matrix for
    the corresponding node numbers.  You pass in the number of the node
    (either as a scalar or an array of locations), and the degrees of
    freedom per node and it returns the corresponding indices in
    the global matrices
    '''
    
    idx = np.array([])
    
    for i in range(len(node)):

        start = DOF*(node[i]-1) + 1
        finish = DOF*node[i]
        stof = np.arange(start,finish+1)
        
        idx = np.append(idx,stof)
    
    # reduce indices by 1 to match pythons order    
    idx = idx - 1
    
    # convert to int
    idx = idx.astype(int)
    
    return idx

def truss(A):
    
    '''
    Compute mass and stress for the 10-bar truss structure
    
    Parameters
    ----------
    A : array of length 10 (row or column)
    cross-sectional areas of each bar.
    See image in book for number order if needed.

    Outputs
    -------
    mass : float
        mass of the entire structure
    stress : array of length 10
        corresponding stress of each bar
 
    '''
    
    P = 1e5 # applied loads
    Ls = 360 # length of sides
    Ld = np.sqrt(2*360**2) # length of diagonals
    
    start = np.array([5, 3, 6, 4, 4, 2, 5, 6, 3, 4])
    finish = np.array([3, 1, 4, 2, 3, 1, 4, 3, 2, 1])
    phi = np.array([0, 0, 0, 0, 90, 90, -45, 45, -45, 45])*np.pi/180
    L = np.array([Ls, Ls, Ls, Ls, Ls, Ls, Ld, Ld, Ld, Ld])
    
    nbar = len(A)
    E = 1e7*np.ones((nbar,)) # modulus of elasticity
    rho = 0.1*np.ones((nbar,)) # material density
    
    Fx = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    Fy = np.array([0.0, -P, 0.0, -P, 0.0, 0.0])
    rigidx = np.array([0, 0, 0, 0, 1, 1])
    rigidy = np.array([0, 0, 0, 0, 1, 1])
    
    n = len(Fx) # number of nodes
    DOF = 2 # number of degrees of freedom
    nbar = len(A) # number of bars
    
    # compute mass
    mass = np.sum(rho*A*L)
    
    # assemble global matrices
    K = np.zeros((DOF*n,DOF*n))
    S = np.zeros((nbar,DOF*n))
    
    for i in range(nbar):
        
        # compute submatrices for the element
        Ksub,Ssub = bar(E[i],A[i],L[i],phi[i])

        # inset submatrix into global matrix
        idx = node2idx([start[i],finish[i]],DOF)
        n_idx = len(idx)
        
        Kind = []
        for ind in range(n_idx):
            Kind.append(idx)
            
        idx_ = np.vstack(Kind)
        K[idx_,idx_.T] = K[idx_,idx_.T] + Ksub
        S[i,idx] = Ssub 

    # setup applied loads
    F = np.zeros((n*DOF,))

    for i in range(n):
        idx = node2idx([i+1],DOF)
        F[idx[0]] = Fx[i]
        F[idx[1]] = Fy[i]
    node_ = np.arange(1,7)    
    idxx = (rigidx != 0);idxx = node_[idxx]
    idxy = (rigidy != 0);idxy = node_[idxy]
    removex = node2idx(idxx,DOF)
    removey = node2idx(idxy,DOF)
    
    nx = len(removex)
    ind_nx = np.arange(0,nx,2)
    removex = removex[ind_nx]
    
    ny = len(removey)
    ind_ny = np.arange(1,ny,2)
    removey = removey[ind_ny]
    
    remove = np.hstack([removex,removey])
    
    K = np.delete(K,remove,axis = 0)
    K = np.delete(K,remove,axis = 1)
    F = np.delete(F,remove)
    S = np.delete(S,remove,axis = 1)
    
    # solve for deflections
    d = solve(K,F)
    
    # compute stress
    stress = S.dot(d)
    
    return mass,stress
    
    
 

