#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def gauss_inter_solve(A, b, x0= None, tol=1e-8, alg='seidel', max_it = 1000):
    """
    

    Parameters
    ----------
    A : array-like, shape = (n,n)
        This contains the coefficient matrix 
    b : array-like, shape = (n,m) with m ≥1 
        This is the right-hand-side vector(s)
    x0 : (optional) array-like, shape = (n,m) or shape = (n,)
        This is an optional kwarg with a default value of None. If no x0 value is imputted 
        the inital guess is taken as all zeros and results in an array. 
        
        If x0 is provided it can either be the same shape as b or a single column vector
        with the n rows. 
        If it is the same shape as b: use x0 directly to find the ouput array.
        If it is a column vector: repeat to initialise an output array with m columns.
    
    This contains the initial guess(es). 
    
    tol : (optional) float
        Relative error tolerance (stopping criterion)
    alg : (optional) string
        This is a flag for the algorithm to be used

    Returns
    -------
    np.ndarray, shape (n,m)
    This is the solution to the system array 
    
    Raises 
    -------
    Value Error: 
    This is the error that will be shown if alg contains a string that does not say which 
    method to use. 
    This is also shown if the shape of A, b and/or x0 do not meet the criterion set.
    
    Type Error:
    This is the error shown if string-like methods are used on alg. 
    
    RuntimeWarning: 
    If the system does not converge before a previously set number of iterations this 
    will alert the user and stop the algorithm. 
    
    -------

    """
    #check A and b are both array-like
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    # check that A is 2d
 
    if (ndimA := len(A.shape)) != 2:
        raise ValueError(f"A is {ndimA}-dimensional, should be 2d")
    # check that A is square
    if not (n := A.shape[0]) == (m := A.shape[1]):
        raise ValueError(f"A has {n} rows and {m} columns, should be square")
        
    # check that b is 1d or 2d
    if (ndimb := len(b.shape)) not in [1, 2]:
        raise ValueError(f"b is {ndimb}-dimensional, should be 1d or 2d")
    nb = b.shape[0]
    if rhs_1d := (ndimb == 1): b = np.reshape(b, (nb, mb := 1))
    else: mb = b.shape[1]    
    # check that b has same number of rows as A
    if (mb := b.shape[0]) != n:
        raise ValueError(f"A has {n} rows and b has {mb} rows, "
                         + "should be equal")                  
    #check alg is either seidel or jacobi
    alg = alg.strip().lower()
   
    if alg not in ['seidel', 'jacobi']:
        raise ValueError ('alg contains a string other than seidel or jacobi')
   
    #set the maximum number of iterations to go to before stopping the calculation
    max_it = 100  
   
    #defining x0
    #make it have the same number of rows as A
    n = len(A[0])
    #make it have the same number of columns as b - if there are more than one
   # m=len(b[1])
    if x0==None:
       x0=np.zeros(n).reshape(-1, 1) #create a column vecor with the same number of rows as A  
               # ensure that x0 is 2d
    if len(x0.shape) == 1: x0 = np.reshape(x0, (n, 1))
        # try to broadcast single initial guess vector to shape of b
    if not rhs_1d and x0.shape[1] == 1:
       x0 = x0 @ np.ones((1, mb))
 
   
        #create a matrix with only the main diagonal
    Ad = np.diag(np.diag(A))
        #calculate the inverse of the diagonal matirx
    Ad_inv = np.diag(1/(np.diag(A)))
        #create a matrix where the main diagonal is 0 so when iterating through you dont
        #have to loop through i and j separately
    A0 = A - Ad
       
        #calculate the normalised B vector
    b_star = Ad_inv @ b
        #calculate the normalised A0 vector
    A0_star = Ad_inv @ A0
       
        #start the iteration counter
    N=1
       
        #set eps_a = 1
    eps_a = 1
       
        #check if the algorithm to be used is seidel or jacobi
        #this is for when the gauss-seidel method is used
    if alg  =='seidel':
        while np.max(eps_a) > tol and N < max_it:
                x_old = np.array(x0)
            #iterating through all the different values in the matrix
                for i, _ in enumerate(A):
                   
                    x0[i, :] = b_star[i:i+1, :] - A0_star[i:i+1, :] @ x0    
                dX = x0 - x_old
                #calculate the approximate error
                eps_a=np.linalg.norm(dX, axis=0)/np.linalg.norm(x0, axis=0)
               
   
               
       
        #this is for when the jacobi method is used    
    elif alg == 'jacobi':
        while np.max(eps_a) > tol and N < max_it:
            x_old = np.array(x0)
            x0 = b_star - A0_star@x0
            dX = x0 - x_old
               
                    #calculate the approximate error
            eps_a=np.linalg.norm(dX, axis=0)/np.linalg.norm(x0, axis=0)
           
            #iterating through all different values in the matrix
    x0=np.reshape(x0,(n,n))        
    return x0
    
    
    