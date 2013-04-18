"""
HW2b
Written by: Daniel McInally  AMATH 583 - HW 2
Derived from demo1.py (Leveque 2013)

This program will determine the polynomial that
passes through the three data points
(-1,0) (1,4) and (2,3)
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

def quad_interp(xi,yi):
        """
        Quadratic interpolation.  Compute the coefficients of the polynomial
        interpolating the points (xi[i],yi[i]) for i = 0,1,2.
        Returns c, an array containing the coefficients of
        p(x) = c[0] + c[1]*x + c[2]*x**2.

        """

        # Error handling of input parameters

        error_message = "xi and yi should have type numpy.ndarray"
        assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

        error_message = "xi and yi should have length 3"
        assert len(xi)==3 and len(yi)==3, error_message


        # Set up the linear system

        A = np.vstack([np.ones(3), xi, xi**2]).T
        b = yi


        #solving the system
        c = solve(A,b)

        return c


def plot_quad (xi,yi):
        """
        Imputs two arrays length 3, computes the coefficient
        matrix from the interpolation, and saves a plot file
        quadratic.png
        """

        c = quad_interp(xi,yi)

        # Plot the polynomial
        x = np.linspace(xi.min() - 1,  xi.max() + 1, 1000)
        y = c[0] + c[1]*x + c[2]*x**2
        
        plt.figure(1)
        plt.clf()
        plt.plot(x,y,'b-')
        
        plt.plot(xi,yi, 'ro')
        plt.ylim(-2,8)
        
        plt.title("Data points and interpolating polynomial")

        plt.savefig('quadratic.png')

def cubic_interp(xi,yi):
        """
        Quadratic interpolation.  Compute the coefficients of the polynomial
        interpolating the points (xi[i],yi[i]) for i = 0,1,2,3.
        Returns c, an array containing the coefficients of
        p(x) = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3

        """

        # Error handling of input parameters

        error_message = "xi and yi should have type numpy.ndarray"
        assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

        error_message = "xi and yi should have length 4"
        assert len(xi)==4 and len(yi)==4, error_message


        # Set up the linear system

        A = np.vstack([np.ones(4), xi, xi**2, xi**3]).T
        b = yi


        #solving the system
        c = solve(A,b)

        return c

def plot_cubic (xi,yi):
        """
        Imputs two arrays length 3, computes the coefficient
        matrix from the interpolation, and saves a plot file
        cubic.png
        """

        c = cubic_interp(xi,yi)

        # Plot the polynomial
        x = np.linspace(xi.min() - 1,  xi.max() + 1, 1000)
        y = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3
        
        plt.figure(1)
        plt.clf()
        plt.plot(x,y,'b-')
        
        plt.plot(xi,yi, 'ro')
        plt.ylim(yi.min()-1,yi.max()+1)
        
        plt.title("Data points and interpolating polynomial")

        plt.savefig('cubic.png')

def poly_interp(xi,yi):
        """
        polynomial interpolation.  Compute the coefficients of the polynomial
        interpolating the points (xi[i],yi[i]) for i = 0,1,2,3,...n.
        Returns c, an array containing the coefficients of
        p(x) = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3 + ... + c[n]*x**n

        """

        # Error handling of input parameters

        error_message = "xi and yi should have type numpy.ndarray"
        assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

        error_message = "xi and yi should have equal length"
        assert len(xi)==len(yi), error_message


        # Set up the linear system
        n=len(xi)
        A=np.array([np.ones(n)], dtype=float)
        for i in range(1,n):
                Ai=np.array(xi**i, dtype=float)
                A = np.vstack([A,Ai])
        A=A.T
        b = yi
        

        #solving the system
        c = solve(A,b)

        return c

def plot_poly (xi,yi):
        """
        Imputs two arrays of any length, and computes the coefficient
        matrix from the interpolation, and saves a plot file
        poly.png
        """

        c = poly_interp(xi,yi)
        n=len(xi)

        # Plot the polynomial
        x = np.linspace(xi.min() - 1,  xi.max() + 1, 1000)
        y = c[n-1]
        for j in range(n-1, 0, -1):
            y = y*x + c[j-1]

        plt.figure(1)
        plt.clf()
        plt.plot(x,y,'b-')
        
        plt.plot(xi,yi, 'ro')
        plt.ylim(yi.min()-1,yi.max()+1)
        
        plt.title("Data points and interpolating polynomial")

        plt.savefig('poly.png')
        
        
# ------------------------------------------------

def test_quad1():
        """
        Test code, no return value or exception if test runs properly.
        """
        xi = np.array([-1.,  0.,  2.])
        yi = np.array([ 1., -1.,  7.])
        c = quad_interp(xi,yi)
        c_true = np.array([-1.,  0.,  2.])
        print "c =      ", c
        print "c_true = ", c_true
        # test that all elements have small error:
        assert np.allclose(c, c_true), \
                "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_quad2():
        """
        Test code, no return value or exception if test runs properly.
        """
        xi = np.array([-5.,  1.,  4.])
        yi = np.array([ 15., 3.,  24.])
        c = quad_interp(xi,yi)
        c_true = np.array([0.,  2.,  1.])
        print "c =      ", c
        print "c_true = ", c_true
        # test that all elements have small error:
        assert np.allclose(c, c_true), \
                "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_cubic1():
        """
        Test code, no return value or exception if test runs properly.
        """
        xi = np.array([ 0.,  1.,  -1., -3.])
        yi = np.array([ 5., 9., 5., -7.])
        c = cubic_interp(xi,yi)
        c_true = np.array([5.,  1.,  2., 1.])
        print "c =      ", c
        print "c_true = ", c_true
        # test that all elements have small error:
        assert np.allclose(c, c_true), \
                "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_poly1():
        """
        Test code, no return value or exception if test runs properly.
        """
        xi = np.array([ 0.,  1.,  -1., -3.])
        yi = np.array([ 5., 9., 5., -7.])
        c = poly_interp(xi,yi)
        c_true = np.array([5.,  1.,  2., 1.])
        print "c =      ", c
        print "c_true = ", c_true
        # test that all elements have small error:
        assert np.allclose(c, c_true), \
                "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_poly2():
        """
        Test code, no return value or exception if test runs properly.
        """
        xi = np.array([ -2.,  -1.,  0., 2., 3.])
        yi = np.array([ 40., 8., 2., 20., 80.])
        c = poly_interp(xi,yi)
        c_true = np.array([2.,  -1.,  3., -1., 1.])
        print "c =      ", c
        print "c_true = ", c_true
        # test that all elements have small error:
        assert np.allclose(c, c_true), \
                "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

# -------------------------------------------------

if __name__=="__main__":
    # "main program"

    print "Running test..."
    test_quad1()
        


