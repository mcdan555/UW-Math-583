"""
    newton.py
    Written by: Daniel McInally
    For: AMATH 583 HW 3

    This module contians the function solve(fvals,x0,debug),
    where fvals(x) returns f(x) and f'(x), x0 is the initial guess,
    and debug is optional flag with default false.

    note: fvals_sqrt and test1 initial was copied from instructor website
    
"""

import numpy as np

def solve(fvals,x0,debug=False)
    """
    This function will use Newton's method to find the roots of the
    equation defined by fvals. fvals contains the definition of the equation
    and will return f(x) and f'(x) at point x.
    """

    # set initial guess and number of iterations

    x = x0
    kmax = 100
    tol = 1.e-14

    if debug:
        print "Initial guess: x = %22.15e" % x

    # do the newton iteration

    for k in range(kmax):
        (f,fp)=fvals(x)
        s0=s
        s=s-f/fp



        if debug:
            print "After %k interations, x = %22.15e" % (k,s)
        
            
    







def fvals_sqrt(x):
    """
    Return f(x) and f'(x) for applying Newton to find a square root.
    """
    f = x**2 - 4.
    fp = 2.*x
    return f, fp

def test1(debug_solve=False):
    """
    Test Newton iteration for the square root with different initial
    conditions.
    """
    from numpy import sqrt
    for x0 in [1., 2., 100.]:
        print " "  # blank line
        x,iters = solve(fvals_sqrt, x0, debug=debug_solve)
        print "solve returns x = %22.15e after %i iterations " % (x,iters)
        fx,fpx = fvals_sqrt(x)
        print "the value of f(x) is %22.15e" % fx
        assert abs(x-2.) < 1e-14, "*** Unexpected result: x = %22.15e"  % x
