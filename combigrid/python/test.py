# -*- coding: utf-8 -*-

from pysgpp import *
import math
import matplotlib.pyplot as plt
import numpy as np
import plot as p

def f(x):
    return math.exp(-x*x)

f = np.vectorize(f)

def copyToVector(list, vec):
    for x in list:
        vec.push_back(x)
    return vec

def testFullGrid():
    interpolator = BarycentricInterpolationEvaluator()
    
    xValues = np.linspace(0, 1, 20)
    plotXValues = np.linspace(0, 1, 200)
    
    yValues = f(xValues)
    
    interpolator.setGridPoints(copyToVector(xValues, DoubleVector()))
    
    plotYValues = []
    
    for v in plotXValues:
        interpolator.setParameter(PyFloatScalarVector(v))
        result = interpolator.eval(copyToVector([PyFloatScalarVector(y) for y in yValues], FloatScalarVector()))
        plotYValues.append(float(result.getValue()))
        
    #Figure erstellen
    fig = plt.figure()
    
    #Achsen mit jeweils 10% Rand erstellen
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.set_xlabel('$x$')
    axes.set_ylabel('$y$')
    axes.set_title('Funktion und Interpolant')
    
    #\, erzeugt ein Leerzeichen in LATEX
    axes.plot(plotXValues, np.array(plotYValues), 'r', label='Interpolant')
    axes.plot(plotXValues, f(plotXValues), 'b', label='$f(x)$')
    axes.legend()
    
    plt.show()


def tensorSquare(vec):
    result = 1.0
    for i in range(len(vec)):
        result *= vec[i] * vec[i]
    return result

def testSingleOperation():
    op = CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(2, multiFunc(tensorSquare))
    
    vec = DataVector(2)
    vec[0] = 0.1
    vec[1] = 0.1
    
    for q in range(0, 11):
        print(q, op.evaluate(q, vec))

class LejaPlotter:
    def __init__(self):
        self.pd = LejaPointDistribution()
        self.evaluator = BarycentricInterpolationEvaluator()
    
        self.numPoints = 25
        self.points = [self.pd.compute(self.numPoints, i) for i in range(self.numPoints)]
        self.pointVector = DoubleVector()
        copyToVector(self.points, self.pointVector)
        self.evaluator.setGridPoints(self.pointVector);
    
    def plot(self):
        p.plot2D(self.getFunc(0), n=500)
        p.plot2D(self.getFunc(5), n=500)
        p.plot2D(self.getFunc(10), n=500)
        p.plot2D(self.getFunc(20), n=500)
        
    def getFunc(self, j):
        return lambda x, this=self, j=j: this.eval(x, j)
        
    def eval(self, x, j):
        self.evaluator.setParameter(PyFloatScalarVector(x))
        
        vec = DoubleVector()
        for i in range(self.numPoints):
            vec.push_back(1.0 if i == j else 0.0)
        
        return self.evaluator.eval(vec).getValue()
    
    
plot = LejaPlotter()
plot.plot()
        
        
#testSingleOperation()


