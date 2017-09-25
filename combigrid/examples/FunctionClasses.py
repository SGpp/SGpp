from sympy import *
import pysgpp
import numpy as np


def displace(x, i, h):
    x_temp = pysgpp.DataVector(x)
    x_temp[i] = x[i] + h
    return x_temp


# get the derivatives for point x
def grad_x(x, function):
    h = 1e-09
    gy = np.zeros(x.shape, dtype=np.float64)

    for k in range(x.shape[0]):
        gy[k] = (function(displace(x, k, h)) -
                 function(displace(x, k, -h))) / (2 * h)

    return gy


"""
x evaluationpoint

returns the k-th gradient with finite differences method
"""


def grad_xk(x, k, function):
    h = 1e-06

    if x[k] <= h:
        gy = (function(displace(x, k, h)) - function(x)) / (h)

    elif x[k] >= 1 - h:

        gy = (function(x) - function(displace(x, k, -h))) / (h)
    else:

        gy = (function(displace(x, k, h)) - function(displace(x, k, -h))) / (2 * h)
    return gy


def getfuncwrapper(func):
    def function(x):
        return func.evalUndisplaced(x)

    return function


# get gradfunction in shape f(x)
def getgradkfunc(func, k):
    def gradk(x):
        return grad_xk(x, k, func)

    return gradk




class funcGradientCollection:
    def __init__(self, func, dim):
        self.func = func
        self.dim = dim

    def getFunction(self):
        return self.func

    def getGradient(self, grad_index_list):
        grad_temp = self.func
        for x in grad_index_list:
            grad_temp = getgradkfunc(grad_temp, x)

        return grad_temp


# includes the function and can create all gradients and mixed gradients with sympy
# naming convention: all variables have to be named x0, x1,...,xn
# TODO maybe get the list of used symbols as alternative

class funcGradientCollectionSymbolic:
    def __init__(self, expr, dim, symbolicDiff=false):
        self.dim = dim
        self.expr = expr

        # name all variables x0,x1....
        self.symbolic_variables = self.init_symbolic_variables()

    def init_symbolic_variables(self):
        symbolic_variables = []
        for i in range(self.dim):
            symbolic_variables.append(Symbol('x' + str(i)))

        return symbolic_variables

    def getGradientSymbolic(self, grad_index_list):

        symbolic_variables = []

        for x in grad_index_list:
            symbolic_variables.append(self.symbolic_variables[x])

        return diff(self.expr, *symbolic_variables)

    def getFunctionSymbolic(self):
        return self.expr
    # return a function that substitutes all symbols with the x value and evaluates
    def getFunction(self):
        symbolic_variables = self.symbolic_variables
        expr = self.expr
        function= lambdify(symbolic_variables, expr, "numpy")

        def func(x):

            x_list=[]
            for i in range(len (x)):
                x_list.append(x[i])

            return function(*x_list)


        return func

    # return a function that substitutes all symbols with the x value and evaluates
    # grad_index_list: list with directions to derivate
    def getGradient(self, grad_index_list):
        symbolic_variables = self.symbolic_variables
        grad_symbolic = self.getGradientSymbolic(grad_index_list)
        gradient = lambdify(symbolic_variables, grad_symbolic, "numpy")
        def grad(x):

          x_list=[]
          for i in range(len(x)):
              x_list.append(x[i])

          return gradient(*x_list)


        return grad


def f(x):
    return x[0] ** 3 * x[1]


#x0, x1, x2 = symbols("x0 x1 x2")
#expr = x0 ** 3 * x1

#test = funcGradientCollectionSymbolic(expr, 2)

#test2 = funcGradientCollection(f, 2)

#print test2.getGradient([0, 1])([1, 1])
#print(test.getGradientSymbolic([0, 0, 1]))
#print(test.getGradient([0, 1])([1, 1]))
