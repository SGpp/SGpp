# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import pysgpp
import math

class TitleFunction(pysgpp.OptObjectiveFunction):
    """Example objective function from the title of my Master's thesis."""
    def __init__(self):
        super(TitleFunction, self).__init__(2)
    
    def eval(self, x):
        """Evaluates the function."""
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])

class TitleFunctionGradient(pysgpp.OptObjectiveGradient):
    """Gradient of TitleFunction."""
    def __init__(self):
        super(TitleFunctionGradient, self).__init__(2)
    
    def evalGradient(self, x, gradient):
        """Evaluates the function gradient."""
        gradient[0] = 8.0 * math.cos(8.0 * x[0])
        gradient[1] = 7.0 * math.cos(7.0 * x[1])
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])

class TitleFunctionHessian(pysgpp.OptObjectiveHessian):
    """Gradient/Hessian of TitleFunction."""
    def __init__(self):
        super(TitleFunctionHessian, self).__init__(2)
    
    def evalHessian(self, x, gradient, hessian):
        """Evaluates the function Hessian."""
        gradient[0] = 8.0 * math.cos(8.0 * x[0])
        gradient[1] = 7.0 * math.cos(7.0 * x[1])
        hessian.set(0, 0, -64.0 * math.sin(8.0 * x[0]))
        hessian.set(0, 1, 0.0)
        hessian.set(1, 0, 0.0)
        hessian.set(1, 1, -49.0 * math.sin(7.0 * x[1]))
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])

class SphereFunctionGradient(pysgpp.OptObjectiveGradient):
    def __init__(self, d):
        super(SphereFunctionGradient, self).__init__(d)
    
    def evalGradient(self, x, gradient):
        """Evaluates the function Gradient."""
        d = self.getDimension()
        fx = 0.0
        for t in range(d):
            y = 10 * x[t] - 1
            fx += y*y
            gradient[t] = 20*y
        return fx

class SphereFunctionHessian(pysgpp.OptObjectiveHessian):
    def __init__(self, d):
        super(SphereFunctionHessian, self).__init__(d)
    
    def evalHessian(self, x, gradient, hessian):
        """Evaluates the function Hessian."""
        d = self.getDimension()
        fx = 0.0
        for t in range(d):
            y = 10 * x[t] - 1
            fx += y*y
            gradient[t] = 20*y
            hessian.set(t, t, 200)
        return fx
