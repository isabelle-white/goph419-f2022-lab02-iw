"""
This script is a test to ensure gauss_inter_solve converges on the correct result when compared to the numpy function 
linalg.solve when given a column vector as the right hand side vector 

Parameters: 
    A: 3x3 matrix used in the gauss iteration 
    b: column vector used as the right hand side vector 
    x0: solution vector produced by gauss_inter_solve 

@author: izzywhite
"""
from linalg_interp_new1 import gauss_inter_solve
import numpy as np

#set up an A matrix which is diagonally dominant for the test to garuantee convergence
A = np.array([[5, -1, -1], 
              [-8, 14, -2],
              [-1, 2, 6]])

#set up the right hand side matrix 
b = np.array([[4.6],
             [15.2],
             [5.0]])

#check that A and b are the right shape and size 
print(A)
print(b)

#solve the linear system using gauss_inter_solve 
x0 = gauss_inter_solve(A, b, x0= None, tol=1e-8, alg='seidel')

#print the solution matrix 
print('The solution is = \n', x0, '\n using the gauss_inter_solve function')
print('''-------''')
#print the solution using the built in numpy function 
print('The solution is = \n', np.linalg.solve(A, b), '\n using the built in numpy function')



"""
This script is a test to ensure gauss_inter_solve produces results which agree with the matrices operation 
A*A^-1 = I 
This test has been completed using both seidel or jacobi iterative methods. 


Parameters: 
    A: 3x3 matrix used in the gauss iteration 
    b: identity matrix of a 3x3 matrix 
    x0_j: solution vector produced by gauss_inter_solve using the jacobi method 
    x0_s: solution vector produced by gauss_inter_solve using the seidel method

@author: izzywhite
"""


from linalg_interp_new import gauss_inter_solve_new
import numpy as np

#set up an A matrix which is diagonally dominant for the test to garuantee convergence
A = np.array([[5, -1, -1],
              [-8, 14, -2],
              [-1, 2, 6]])
 
#set up the right hand side matrix
b = np.eye(3)
 
#check that A and b are the right shape and size
print(A)
print(b)
print('----')
 
#solve the linear system using gauss_inter_solve
x0_j = gauss_inter_solve_new(A, b, x0= None, tol=1e-8, alg='jacobi')
 
#print the solution matrix
print('The solution is = \n', x0_j, '\n using the gauss_inter_solve function - Jacobi')
#print the solution using the built in numpy function
print('The solution is = \n', np.linalg.solve(A, b), '\n using the built in numpy function')
#check A * A inverse (which is x0) = the identity
rhs = A@x0_j
print('A @ x0 = \n', rhs)

print('----')

#solve the linear system using gauss_inter_solve
x0_s = gauss_inter_solve_new(A, b, x0= None, tol=1e-8, alg='seidel')
 
#print the solution matrix
print('The solution is = \n', x0_s, '\n using the gauss_inter_solve function - Seidel')
#print the solution using the built in numpy function
print('The solution is = \n', np.linalg.solve(A, b), '\n using the built in numpy function ')
#check A * A inverse (which is x0) = the identity
rhs = A@x0_s
print('A @ x0 = \n', rhs)


"""
This script generates data from linear, quadratic, and cubic functions using the spline function. 
On the same plot it compares the data produced by the origional function to ensure the spline function is working 
fully. 

Parameters: 
    xd: array-like like of 10 evenly spaced floats from -10 to 10 
    y_liner = linear function using xd
    y_quadratic = quadratic function using xd
    y_cubic = cubic function using xd
    f_1 = 1st order spline function 
    y_1 = y values using order = 1 spline
    f_2 = 2nd order spline function
    y_2 = y values using order = 2 spline
    f_3 = 3rd order spline function 
    y_2 = y values using order = 3 spline

Inputs:
    A linear, quadratic, and cubic function 
    
Outputs: 
    Plot with each spline function and respective true data plotted in a different colour

@author: izzywhite
"""
import numpy as np
from spline_function_new import spline_function
import matplotlib.pyplot as plt 

#create an array of 10 x values from 1-10
xd = np.linspace(-10, 10, 10)

#create linear, quadratic and cubic functions to use for spline function 
y_linear = 50 *xd

y_quadratic = 5*(xd**2) + 3*xd + xd+2

y_cubic = xd**3 + 3*(xd**2) - 6*xd - (xd-8)



#calculate the y values using the spline function and x data inputted  
f_1 = spline_function(xd, y_linear, order=1)
y_1 = f_1(xd) 

f_2 = spline_function(xd, y_quadratic, order=2)
y_2 = f_2(xd)

f_3 = spline_function(xd, y_cubic, order=3)
y_3 = f_3(xd)

#plot the different functions on the same graph adding markers for the values calculated using the exact function 
plt.plot(xd, y_linear, 'gD', markersize = '4')
plt.plot(xd, y_1, 'g', label ="50*xd")

plt.plot(xd, y_quadratic, 'ro', markersize = '4')
plt.plot(xd, y_2, 'r', label = 'xd**2 + 3*xd + 2')

plt.plot(xd, y_cubic, 'mx', markersize = '4')
plt.plot(xd, y_3, 'm', label = 'xd**3 + 3*(xd**2) - 6*xd - 8')

#add graph title and axis labels
plt.title("Spline Function for different order systems")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()

"""
This script compares my spline function with the Univariate spline function for a third order polynomial
by taking an expenential function and plotting the results of each on 2 separate graphs.

Parameters:
    xd: array-like, evenly spaced floats 
    yd: array-like, data produced by the exponential function 
    f_univ = spline function using the univariate spline
    y_univ: array-like, y values produced by the univariate spline function 
    f_interp: spline function using  spline_function and order=3
    y_new: array-like, y values produced by my spline function 

Input: 
    Exponential function 

Output:
    Graphs of the spline function and Univariate spline function

@author: izzywhite
"""

import numpy as np
from spline_function_new import spline_function
from scipy.interpolate import UnivariateSpline  
import matplotlib.pyplot as plt 

xd = np.linspace(1, 10, 50)
yd = np.exp(xd)

#plot the result when using the built in numpy univariate spline function 
f_univ = UnivariateSpline(xd, yd, k=3, s=0, ext = 'raise')
y_univ = f_univ(xd)
plt.plot(xd, yd, 'ro', markersize ="4", label='raw data')
plt.plot(xd, y_univ, 'g', label='univariate spline')
plt.title("Univariate Spline Function")
plt.legend()
plt.show()

#plot the result when using spline_function 
f_interp = spline_function(xd, yd, order=3)
y_new = f_interp(xd) 

plt.plot(xd, yd, 'ro', markersize ="4", label='raw data')
plt.plot(xd, y_new, 'g', label='spline function')
plt.title("My Spline Function")
plt.legend()
plt.show()

