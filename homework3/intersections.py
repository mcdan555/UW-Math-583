"""

    intersections.py
    Written by: Daniel McInally
    For: AMATH 583 HW 3

    This program will use newton.solve to find the intersections
    of g1(x)=x cos(Pi*x) and g2(x) = 1 - 0.6 * x^2
"""

from newton import solve
from numpy import *
from pyplot import *
from math import *


def fval_g1(x):
    y=x*cos(pi*x)
    return y

def fval_g2 (x):
    y=1-0.6*x**2
    return y

def fval_f (x):
    f=x*cos(pi*x)-1-0.6*x**2
    fp=cos(pi*x)-pi*x*sin(pi*x)-1.2*x
    return f,fp

def plotfuncs (sol):

    x=linspace(-10,10,2000)

    for k in x:
        plot(k,fval_g1(k))
        plot(k,fval_g2(k))

    for k in sol:
        plot(k,fval_g1(k),'b.')

def findintersections:
    """
    This will use newton.solve to find the intersections and return
    the array sol
    """

    x0 = [-2.18026,-1.61399,0.794267,1.44331] #tupple of initial values
    sol=[]

    for k in xo:
        print " "
        x,iters = solve(fvals_f,k,debug=false)
        print "solve returns x = %22.15e after %i interations " % (x,iters)
        f,fp = fvals_f(x)
        print "the value of f(x) is %22.15e" % f
        assert abs(x-2.) < 1e-14, "*** Unexpected result: x = %22.15e" % x
        sol.append[x]
    
    

# --------- Main Program -----------

if __name__=="__main__":

    inter_sol = findintersections()

    plotfuncs (sol)
    

    
    
