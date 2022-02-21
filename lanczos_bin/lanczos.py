import numpy as np
import scipy as sp

def exact_lanczos(A,q0,k,B=None,reorth=True):
    """
    run Lanczos with reorthogonalization
    
    Input
    -----
    A : entries of diagonal matrix A
    q0 : starting vector
    k : number of iterations
    B : entries of diagonal weights for orthogonalization
    """
    
    n = len(A)
    
    if B is None:
        B = np.ones(n,dtype=A.dtype)
    
    Q = np.zeros((n,k),dtype=A.dtype)
    a = np.zeros(k,dtype=A.dtype)
    b = np.zeros(k,dtype=A.dtype)
    
    Q[:,0] = q0 / np.sqrt(q0*B@q0)
    
    for i in range(1,k+1):
        # expand Krylov space

      #  if i>1:
       #     print(b[i-2],qi@Q[:,i-2])

        qi = A*Q[:,i-1] - b[i-2]*Q[:,i-2] if i>1 else A*Q[:,i-1]
        
        a[i-1] = (qi*B)@Q[:,i-1]
        qi -= a[i-1]*Q[:,i-1]
        
        if reorth:
            qi -= Q[:,:i-2]@(Q[:,:i-2].T@(B*qi))
            
        b[i-1] = np.sqrt((qi*B)@qi)
        if i < k:
            Q[:,i] = qi / b[i-1]
                
    return Q,(a,b)

def lanczos_FA(f,Q,a_,b_,normb=1):
    
    try:
        theta,S = sp.linalg.eigh_tridiagonal(a_,b_,tol=1e-30)
    except:
        T = np.diag(a_) + np.diag(b_,-1) + np.diag(b_,1)
        theta,S = sp.linalg.eigh(T)
        
#    fT = S*f(theta)@S.T 
    
    # e0 = Q^Tb in exaxt arithmetic. IDK if this matters numerically..
    e0 = np.zeros_like(Q[0])
    e0[0] = normb
    
    return Q@(S@(f(theta)*(S.T@e0)))

