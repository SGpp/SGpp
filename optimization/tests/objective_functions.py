# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import pysgpp
import math

class ExampleFunction(pysgpp.OptObjectiveFunction):
    """Example objective function from the title of my Master's thesis."""
    def __init__(self):
        super(ExampleFunction, self).__init__(2)
    
    def eval(self, x):
        """Evaluates the function."""
        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
            return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])
        else:
            return float("inf")

class ExampleFunctionGradient(pysgpp.OptObjectiveGradient):
    """Gradient of ExampleFunction."""
    def __init__(self):
        super(ExampleFunctionGradient, self).__init__(2)
    
    def eval(self, x, gradient):
        """Evaluates the function gradient."""
        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
            gradient[0] = 8.0 * math.cos(8.0 * x[0])
            gradient[1] = 7.0 * math.cos(7.0 * x[1])
            return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])
        else:
            return float("inf")

class ExampleFunctionHessian(pysgpp.OptObjectiveHessian):
    """Gradient/Hessian of ExampleFunction."""
    def __init__(self):
        super(ExampleFunctionHessian, self).__init__(2)
    
    def eval(self, x, gradient, hessian):
        """Evaluates the function Hessian."""
        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
            gradient[0] = 8.0 * math.cos(8.0 * x[0])
            gradient[1] = 7.0 * math.cos(7.0 * x[1])
            hessian.set(0, 0, -64.0 * math.sin(8.0 * x[0]))
            hessian.set(0, 1, 0.0)
            hessian.set(1, 0, 0.0)
            hessian.set(1, 1, -49.0 * math.sin(7.0 * x[1]))
            return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])
        else:
            return float("inf")



class SphereFunctionGradient(pysgpp.OptObjectiveGradient):
    def __init__(self, d):
        super(SphereFunctionGradient, self).__init__(d)
    
    def eval(self, x, gradient):
        """Evaluates the function gradient."""
        d = self.getDimension()
        
        if all([0.0 <= x[t] <= 1.0 for t in range(d)]):
            fx = 0.0
            for t in range(d):
                y = 10 * x[t] - 1
                fx += y*y
                gradient[t] = 20*y
            return fx
        else:
            return float("inf")

class SphereFunctionHessian(pysgpp.OptObjectiveHessian):
    def __init__(self, d):
        super(SphereFunctionHessian, self).__init__(d)
    
    def eval(self, x, gradient, hessian):
        """Evaluates the function Hessian."""
        d = self.getDimension()
        
        if all([0.0 <= x[t] <= 1.0 for t in range(d)]):
            fx = 0.0
            for t in range(d):
                y = 10 * x[t] - 1
                fx += y*y
                gradient[t] = 20*y
                hessian.set(t, t, 200)
            return fx
        else:
            return float("inf")



class G3ObjectiveFunction(pysgpp.OptObjectiveFunction):
    def __init__(self, d):
        super(G3ObjectiveFunction, self).__init__(d)
    
    def eval(self, x):
        """Evaluates the function."""
        d = self.getDimension()
        
        if all([0.0 <= x[t] <= 1.0 for t in range(d)]):
            fx = d ** (d / 2.0)
            for t in range(d):
                fx *= x[t]
            return fx
        else:
            return float("inf")

class G3ObjectiveGradient(pysgpp.OptObjectiveGradient):
    def __init__(self, d):
        super(G3ObjectiveGradient, self).__init__(d)
    
    def eval(self, x, gradient):
        """Evaluates the function gradient."""
        d = self.getDimension()
        
        if all([0.0 <= x[t] <= 1.0 for t in range(d)]):
            fx = d ** (d / 2.0)
            
            for t in range(d):
                gradient[t] = fx
            
            for t in range(d):
                for t2 in range(d):
                    if t2 != t:
                        gradient[t2] *= x[t]
                    else:
                        fx = fx * x[t]
            
            return fx
        else:
            return float("inf")



class G3ConstraintFunction(pysgpp.OptConstraintFunction):
    def __init__(self, d):
        super(G3ConstraintFunction, self).__init__(d, 1)
    
    def eval(self, x, value):
        """Evaluates the constraint."""
        d = self.getDimension()
        
        if all([0.0 <= x[t] <= 1.0 for t in range(d)]):
            gx = -1.0
            for t in range(d):
                gx += x[t] * x[t]
            value[0] = gx
        else:
            value[0] = float("inf")

class G3ConstraintGradient(pysgpp.OptConstraintGradient):
    def __init__(self, d):
        super(G3ConstraintGradient, self).__init__(d, 1)
    
    def eval(self, x, value, gradient):
        """Evaluates the constraint gradient."""
        d = self.getDimension()
        
        if all([0.0 <= x[t] <= 1.0 for t in range(d)]):
            gx = -1.0
            for t in range(d):
                gx += x[t] * x[t]
                gradient.set(0, t, 2.0 * x[t])
            value[0] = gx
        else:
            value[0] = float("inf")



#class G6ObjectiveFunction(pysgpp.OptObjectiveFunction):
#    def __init__(self):
#        super(G6ObjectiveFunction, self).__init__(2)
#    
#    def eval(self, x):
#        """Evaluates the function."""
#        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
#            x0 = x[0] * (100.0 - 13.0) + 13.0
#            x1 = x[1] * (100.0 -  0.0) +  0.0
#            fx = (x0 - 10.0)**3 + (x1 - 20.0)**3
#            return fx
#        else:
#            return float("inf")
#
#class G6ObjectiveGradient(pysgpp.OptObjectiveGradient):
#    def __init__(self):
#        super(G6ObjectiveGradient, self).__init__(2)
#    
#    def eval(self, x, gradient):
#        """Evaluates the function gradient."""
#        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
#            x0 = x[0] * (100.0 - 13.0) + 13.0
#            x1 = x[1] * (100.0 -  0.0) +  0.0
#            fx = (x0 - 10.0)**3 + (x1 - 20.0)**3
#            gradient[0] = 3.0 * (x0 - 10.0)**2 * (100.0 - 13.0)
#            gradient[1] = 3.0 * (x1 - 20.0)**2 * (100.0 -  0.0)
#            
#            return fx
#        else:
#            return float("inf")
#
#
#
#class G6ConstraintFunction(pysgpp.OptConstraintFunction):
#    def __init__(self):
#        super(G6ConstraintFunction, self).__init__(2, 2)
#    
#    def eval(self, x, value):
#        """Evaluates the constraint."""
#        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
#            x0 = x[0] * (16.0 - 13.0) + 13.0
#            x1 = x[1] * (12.0 -  0.0) +  0.0
#            value[0] = -(x0 - 5.0)**2 - (x1 - 5.0)**2 + 100.0
#            value[1] =  (x0 - 6.0)**2 + (x1 - 5.0)**2 - 82.81
#        else:
#            value[0] = float("inf")
#            value[1] = float("inf")
#
#class G6ConstraintGradient(pysgpp.OptConstraintGradient):
#    def __init__(self):
#        super(G6ConstraintGradient, self).__init__(2, 2)
#    
#    def eval(self, x, value, gradient):
#        """Evaluates the constraint gradient."""
#        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
#            x0 = x[0] * (16.0 - 13.0) + 13.0
#            x1 = x[1] * (12.0 -  0.0) +  0.0
#            value[0] = -(x0 - 5.0)**2 - (x1 - 5.0)**2 + 100.0
#            value[1] =  (x0 - 6.0)**2 + (x1 - 5.0)**2 - 82.81
#            gradient.set(0, 0, -2.0 * (x0 - 5.0) * (16.0 - 13.0))
#            gradient.set(0, 1, -2.0 * (x1 - 5.0) * (12.0 -  0.0))
#            gradient.set(1, 0,  2.0 * (x0 - 6.0) * (16.0 - 13.0))
#            gradient.set(1, 1,  2.0 * (x1 - 5.0) * (12.0 -  0.0))
#        else:
#            value[0] = float("inf")
#            value[1] = float("inf")



class G8ObjectiveFunction(pysgpp.OptObjectiveFunction):
    def __init__(self):
        super(G8ObjectiveFunction, self).__init__(2)
    
    def eval(self, x):
        """Evaluates the function."""
        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
            x0 = 10.0 * x[0]
            x1 = 10.0 * x[1]
            fx = -math.sin(2.0*math.pi*x0)**3 * math.sin(2.0*math.pi*x1) / \
                 (x0**3 * (x0 + x1))
            return fx
        else:
            return float("inf")

class G8ObjectiveGradient(pysgpp.OptObjectiveGradient):
    def __init__(self):
        super(G8ObjectiveGradient, self).__init__(2)
    
    def eval(self, x, gradient):
        """Evaluates the function gradient."""
        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
            x0 = 10.0 * x[0]
            x1 = 10.0 * x[1]
            fx = -math.sin(2.0*math.pi*x0)**3 * math.sin(2.0*math.pi*x1) / \
                 (x0**3 * (x0 + x1))
            gradient[0] = 6.0 * math.pi * math.cos(2.0*math.pi*x0) * \
                          math.sin(2.0*math.pi*x0)**2 * \
                          math.sin(2.0*math.pi*x1) / ((x0 + x1) * x0**3) - \
                          3.0 * math.sin(2.0*math.pi*x0)**3 * \
                          math.sin(2.0*math.pi*x1) / ((x0 + x1) * x0**4) - \
                          math.sin(2.0*math.pi*x0)**3 * \
                          math.sin(2.0*math.pi*x1) / ((x0 + x1)**2 * x0**3)
            gradient[1] = 2.0 * math.pi * math.cos(2.0*math.pi*x1) * \
                          math.sin(2.0*math.pi*x0)**3 / ((x0 + x1) * x0**3) - \
                          math.sin(2.0*math.pi*x0)**3 * \
                          math.sin(2.0*math.pi*x1) / ((x0 + x1)**2 * x0**3)
            gradient[0] *= -10.0
            gradient[1] *= -10.0
            
            return fx
        else:
            return float("inf")



class G8ConstraintFunction(pysgpp.OptConstraintFunction):
    def __init__(self):
        super(G8ConstraintFunction, self).__init__(2, 2)
    
    def eval(self, x, value):
        """Evaluates the constraint."""
        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
            x0 = 10.0 * x[0]
            x1 = 10.0 * x[1]
            value[0] = x0**2 - x1 + 1.0
            value[1] = 1.0 - x0 + (x1 - 4.0)**2
        else:
            value[0] = float("inf")
            value[1] = float("inf")

class G8ConstraintGradient(pysgpp.OptConstraintGradient):
    def __init__(self):
        super(G8ConstraintGradient, self).__init__(2, 2)
    
    def eval(self, x, value, gradient):
        """Evaluates the constraint gradient."""
        if all([0.0 <= x[t] <= 1.0 for t in range(2)]):
            x0 = 10.0 * x[0]
            x1 = 10.0 * x[1]
            value[0] = x0**2 - x1 + 1.0
            value[1] = 1.0 - x0 + (x1 - 4.0)**2
            gradient.set(0, 0, 2.0 * x0 * 10.0)
            gradient.set(0, 1, -1.0 * 10.0)
            gradient.set(1, 0, -1.0 * 10.0)
            gradient.set(1, 1, 2.0 * (x1 - 4.0) * 10.0)
        else:
            value[0] = float("inf")
            value[1] = float("inf")
