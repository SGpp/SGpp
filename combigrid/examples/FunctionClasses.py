import sympy as sp
import pysgpp
import numpy as np



from itertools import chain, combinations


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))




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
    def __init__(self, expr, dim):
        self.dim = dim
        self.expr = expr

        # name all variables x0,x1....
        self.symbolic_variables = self.init_symbolic_variables()

    def init_symbolic_variables(self):
        symbolic_variables = []
        for i in range(self.dim):
            symbolic_variables.append(sp.Symbol('x' + str(i)))

        return symbolic_variables

    def getGradientSymbolic(self, grad_index_list):

        symbolic_variables = []

        for x in grad_index_list:
            symbolic_variables.append(self.symbolic_variables[x])

        return sp.diff(self.expr, *symbolic_variables)

    def getFunctionSymbolic(self):
        return self.expr
    # return a function that substitutes all symbols with the x value and evaluates
    def getFunction(self):
        symbolic_variables = self.symbolic_variables
        expr = self.expr
        function= sp.lambdify(symbolic_variables, expr, "numpy")

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
        gradient = sp.lambdify(symbolic_variables, grad_symbolic, "numpy")
        def grad(x):

          x_list=[]
          for i in range(len(x)):
              x_list.append(x[i])

          return gradient(*x_list)


        return grad


def f(x):
    return x[0] ** 3 * x[1]

def BraninSymbolic():
    x0, x1, x2 = sp.symbols("x0 x1 x2")

    x1_tmp = 15.0 * x0 - 5.0;

    x2_tmp = 15.0 * x1

    tmp = x2_tmp - 5.1 * x1_tmp * x1_tmp / (4.0 * sp.pi * sp.pi) + 5.0 * x1_tmp / sp.pi - 6.0;

    return tmp * tmp + 10.0 * (1.0 - 1.0 / (8.0 * sp.pi)) * sp.cos(x1_tmp) + 10.0;

def testfSymbolic_2d():
    x0, x1, x2 = sp.symbols("x0 x1 x2")

    x0_temp=5.0*x0
    x1_temp=5.0*x1
    return x0**2*sp.cos(x1_temp)**2*x1+(x0**2+sp.cos(x1_temp))**2

def testfSymbolic_4d():
    x0, x1, x2,x3 = sp.symbols("x0 x1 x2 x3")

    x0_temp=5.0*x0
    x1_temp=5.0*x1
    x2_temp=4.0*x2
    x3_temp=3.0*x3
    return x0**2*sp.cos(x1_temp)**2*x1*x3_temp**4+(x0**2+sp.cos(x1_temp))**2*x2_temp

def testfSymbolic2_4d():
    x0, x1, x2,x3 = sp.symbols("x0 x1 x2 x3")
    x0_temp = 5.0 * x0
    x1_temp = 5.0 * x1
    x2_temp = 4.0 * x2
    x3_temp = 3.0 * x3

    return x0_temp**2*sp.cos(0.5-x1_temp)*x2_temp**3*sp.exp(x3_temp)


def test2DSymbolic():
    x0, x1, x2 = sp.symbols("x0 x1 x2")
    return sp.cos(x1*2)*sp.sin(0.5-x0*4)**2

def get_symbolictest():
    x0, x1, x2 = sp.symbols("x0 x1 x2")
    return (x0 * x1 + 1) ** 2 * sp.sin(x0)

def easySymbolic():
    x0, x1, x2 = sp.symbols("x0 x1 x2")
    return x0 * x1*sp.sin(x0)*sp.sin(x1)

def testfuncione():
    x0, x1, x2 = sp.symbols("x0 x1 x2")

    erg=x1

    for i in range(5):
        erg+=x0
    return erg


def RosenbrockSymbolic(d):

  result = 0.0;

  x=[]

  for i in range (d):
    x.append(sp.symbols("x"+str(i)))


  xt = 15.0 * x[0] - 5.0;

  for t in range(1,d):
    xtm1 = xt;
    xt = 15.0 * x[t] - 5.0;

    tmp1 = xt - xtm1 * xtm1;
    tmp2 = 1.0 - xtm1;
    result += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;


  return result;

def f1D_test(x):
    return 1


#x0, x1, x2 = symbols("x0 x1 x2")
#expr = x0 ** 3 * x1

#test = funcGradientCollectionSymbolic(expr, 2)

#test2 = funcGradientCollection(f, 2)

#print test2.getGradient([0, 1])([1, 1])
#print(test.getGradientSymbolic([0, 0, 1]))
#print(test.getGradient([0, 1])([1, 1]))
