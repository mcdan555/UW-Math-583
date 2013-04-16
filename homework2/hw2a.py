"""
Written by: Daniel McInally  AMATH 583 - HW 2
Derived from demo1.py (Leveque 2013)

This program will determine the polynomial that
passes through the three data points
(-1,0) (1,4) and (2,3)
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

# Set up the linear system

# Data points
xi = np.array([-1.,1.,2.])
yi = np.array([0.,4.,3.])

# matrix and b value
A = np.array([[1.,-1.,1.],[1.,0.,0.],[1.,2.,4.]])
b = yi


#solving the system
c = solve(A,b)

print "The polynomial coefficients are:"
print c

# Plot the polynomial
x = np.linspace(-2,3,1001)
y = c[0] + c[1]*x + c[2]*x**2

plt.figure(1)
plt.clf()
plt.plot(x,y,'b-')

plt.plot(xi,yi, 'ro')
plt.ylim(-2,8)

plt.title("Data points and interpolating polynomial")

plt.savefig('hw2a.png')


