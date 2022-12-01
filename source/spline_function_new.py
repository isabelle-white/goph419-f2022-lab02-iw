#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def spline_function(xd, yd, order=3):
    """
    Created on Sat Nov 12 15:36:32 2022

    Parameters 
    ----------
    xd: array-like of float 
    yd: array-like of float with the same shape as xd
    order: this is an optional integer of either 1, 2, or 3 with a default of 3
    
    Returns
    ------
    function that takes one parameter and returns an interpolated value y
    
    Raises
    ------
    Value Error:
    Length of xd and yd are different 
    There are repeated values in xd
    xd values are not in increasing order
    order is a value other than 1, 2, or 3

    @author: izzywhite
    """
          
    xd = np.array(xd, dtype=float)
    yd = np.array(yd, dtype=float)  
        
    if len(xd) != len(yd):
        raise ValueError('number of independent and dependent variables must be equal')
    
    #test to check if xd is non-unique 
    xd_unique = np.unique(xd) 
    if len(xd_unique) != len(xd):
        raise ValueError('there are repeated values in xd')
      
        #test to check that xd is sorted
    if (all(xd[i] >= xd[i+1] for i in range(len(xd)-1))):
        raise ValueError('xd is not sorted')
    
    
    #taking the differences between adjacent x values
    dx = np.diff(xd)
    dy = np.diff(yd)
    #first order divided difference
    df1 = dy/dx
     
    
    if order==1:
        def s1(x):
  
            #finding a coefficient 
            a = yd[:-1]
            
            #finding b coefficient 
            b = df1[:-1]
            
            #finding the index which x is between adjacent xd values
            for xk in x:
                
                #finding the index of xd in x                 
                k = np.array([np.nonzero(xd>=xk)[0][0] - 1 for xk in x])
                
                k = np.where(k<0, 0, k)
            
            #defining spline function
                y = a[k-1] + b[k-1]*(x-xd[k-1]) 
                return y
        
        return s1 
                
        
    elif order ==2:
        def s2(x):
            
            #making the right hand side matrix using the dx out of the elif
            N = xd.shape[0]
            rhs = np.zeros(N-1)
            rhs[1:] = np.diff(df1, axis=0)
            #rhs = np.vstack([[[0]], np.diff(df1, axis=0)]) - not needed
            
            
            #finding the length of the input data array 
            N = len(xd)
            #create an empty matrix of length N by N
            A = np.zeros((N-1,N-1))
            
            #finding the vaues for the first and last rows of the A matrix
            A[0,0:2] = [1, -1]
            A[1:, : -1] += np.diag(dx[:-1])
            A[1:, 1:] += np.diag(dx[1:])
            
            #finding the c coefficients by solving the linear system using gauss 
            #elimination function from question 1
            c = np.linalg.solve(A, rhs)
        
        
            #find the b coefficient for the function 
            b = dy - c*dx 
            #find the a coefficient for the function 
            a = yd[:-1]
            
            #finding the index which x is between adjacent xd values
            for xk in x:
                
                #finding the index of xd in x                 
                k = np.array([np.nonzero(xd>=xk)[0][0] - 1 for xk in x])
                
                k = np.where(k<0, 0, k)
            
            #defining spline function
                y = a[k] + b[k]*(x-xd[k]) + c[k]*(x-xd[k])**2 
                
                return y
    

        return s2
                
    
    elif order==3:
       def s3(x):

           #making the right hand side matrix using dx from outside elif
           #rhs = np.vstack([[[0]], 3*np.diff(df1, axis=0), [[0]]])

           #finding the length of the input data array 
           N = xd.shape[0]
           #create an empty matrix of length N by N
           d2 = np.diff(df1)
           rhs = np.zeros(N)
           rhs[1:-1] = 3. * d2
           A = np.zeros((N,N))
           #finding the vaues for the first and last rows of the A matrix
           A[1:-1, 1:-1] = (2. * (np.diag(dx[:-1]) + np.diag(dx[1:]))
                            + np.diag(dx[1:-1], -1) + np.diag(dx[1:-1], 1))
           A[1, 0] = dx[0]
           A[-2, -1] = dx[-1]
           A[0, :3] = [dx[1], -(dx[0] + dx[1]), dx[0]]
           A[-1, -3:] = [dx[-1], -(dx[-2] + dx[-1]), dx[-2]]
           #finding the values for the 'bulk' of the A matrix 
           #finding the c coefficients by solving the linear system using gauss 
           #elimination function from question 1
           c = np.linalg.solve(A,rhs)
           #find the d coefficient for the function 
           d = np.diff(c) / (3. * dx)
           #find the b coefficient for the function 
           b = df1 - dx * (2 * c[:-1] + c[1:]) / 3
           
           
          #finding the index which x is between adjacent xd values
           k = np.array([np.nonzero(xd >= xk)[0][0] - 1 for xk in x])
           k = np.where(k < 0, 0, k)
           y = np.array([(yd[k]
                          + b[k] * (xk - xd[k])
                          + c[k] * (xk - xd[k]) ** 2
                          + d[k] * (xk - xd[k]) ** 3)
                         for k, xk in zip(k, x)])
           return y
       return s3
        
     
    else:
        raise ValueError('order given is not 1, 2, or 3')
        
   
        