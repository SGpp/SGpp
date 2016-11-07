#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page combigrid_gettingStarted_py gettingStarted.py (Start Here)

## At the beginning of the program, we have to import the pysgpp library.
import pysgpp
import math


## Let's define a function for evaluation
def f(x):
    product = 1.0
    for i in range(x.getSize()):
        product *= math.exp(-x[i])
    return product

func = pysgpp.multiFunc(f)
d = 3


def example1():
    growthFactor = 2

    operation = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, func, growthFactor)

    result = operation.evaluate(2)

    print("Quadrature result: " + str(result) + ", analytical solution: " + str(math.pow(1.0 - 1.0/math.e, d)))

    print("Number of function evaluations: " + str(operation.getStorage().getNumEntries()))


def example2():
    operation = pysgpp.CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(d, func)

    evaluationPoint = pysgpp.DataVector([0.1572, 0.6627, 0.2378])

    result = operation.evaluate(3, evaluationPoint)

    print("Interpolation result: " + str(result) + ", function value: " + str(func(evaluationPoint)))

    print("Function evaluations: " + str(operation.getStorage().getNumEntries()))

    evaluationPoint[0] = 0.4444
    print("Target function value: " + str(func(evaluationPoint)))
    operation.setParameters(evaluationPoint)

    levelManager = operation.getLevelManager()

    levelManager.addRegularLevels(3)

    print("Regular result 1: " + str(operation.getResult()))
    print("Total function evaluations: " + str(operation.getStorage().getNumEntries()))

    levelManager.addRegularLevelsByNumPointsParallel(50, 4)
    print("Regular result 2: " + str(operation.getResult()))
    print("Total function evaluations: " + str(operation.getStorage().getNumEntries()))

    operation.setLevelManager(pysgpp.AveragingLevelManager())
    levelManager = operation.getLevelManager()

    levelManager.addLevelsAdaptive(60)
    print("Adaptive result: " + str(operation.getResult()))
    print("Total function evaluations: " + str(operation.getStorage().getNumEntries()))


def example3():
    operation = pysgpp.CombigridMultiOperation.createLinearLejaPolynomialInterpolation(d, func)

    firstParam = [0.2, 0.6, 0.7]
    secondParam = [0.3, 0.9, 1.0]
    parameters = pysgpp.DataMatrix([firstParam, secondParam])

    stopwatch = pysgpp.Stopwatch()
    result = operation.evaluate(3, parameters)
    stopwatch.log()
    print("First result: " + str(result[0]) + ", function value: " + str(func(pysgpp.DataVector(firstParam))))
    print("Second result: " + str(result[1]) + ", function value: " + str(func(pysgpp.DataVector(secondParam))))


def loggingF(x):
    print("call function")
    return x[0]


def example4():
    loggingFunc = pysgpp.multiFunc(loggingF)
    lookupTable = pysgpp.FunctionLookupTable(loggingFunc)

    operation = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, lookupTable.toMultiFunction())

    result = operation.evaluate(2)
    print("Result computed: " + str(result))

    pysgpp.writeToFile("lookupTable.log", lookupTable.serialize())

    restoredLookupTable = pysgpp.FunctionLookupTable(loggingFunc)
    restoredLookupTable.deserialize(pysgpp.readFromFile("lookupTable.log"))
    operation2 = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, restoredLookupTable.toMultiFunction())

    result = operation2.evaluate(2)
    print("Result computed (2nd time): " + str(result))

    pysgpp.writeToFile("storage.log", operation.getStorage().serialize())
    operation3 = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, pysgpp.multiFunc(loggingFunc))
    operation3.getStorage().deserialize(pysgpp.readFromFile("storage.log"))
    result = operation3.evaluate(2)
    print("Result computed (3rd time): " + str(result))


def example5():
    grids = pysgpp.AbstractPointHierarchyVector()
    grids.push_back(pysgpp.CombiHierarchies.expChebyshev())
    grids.push_back(pysgpp.CombiHierarchies.linearLeja(3))
    grids.push_back(pysgpp.CombiHierarchies.expUniform())

    evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
    evaluators.push_back(pysgpp.CombiEvaluators.polynomialInterpolation())
    evaluators.push_back(pysgpp.CombiEvaluators.quadrature())
    evaluators.push_back(pysgpp.CombiEvaluators.linearInterpolation())

    levelManager = pysgpp.WeightedRatioLevelManager()
    operation = pysgpp.CombigridOperation(grids, evaluators, levelManager, func)

    parameters = pysgpp.DataVector([0.777, 0.14159])
    result = operation.evaluate(2, parameters)
    print("Result: " + str(result))


# Call the examples

print("Example 1:")
example1()

print("\nExample 2:")
example2()

print("\nExample 3:")
example3()

print("\nExample 4:")
example4()

print("\nExample 5:")
example5()
