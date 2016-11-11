#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page combigrid_gettingStarted_py gettingStarted.py (Start Here)
## This tutorial contains examples with increasing complexity to introduce you to the combigrid
## module. The combigrid module is quite separated from the other modules. It only refers to the
## base module for things like DataVector and DataMatrix.

## At the beginning of the program, we have to import the pysgpp library.
import pysgpp
import math


## The first thing we need is a function to evaluate. This function will be evaluated on the domain
## \f$[0, 1]^d\f$. This particular function can be used with any number of dimensions.
## The input parameter of the function is of type pysgpp.DataVector, so do not treat it like a list.
## The return type is float.
def f(x):
    product = 1.0
    for i in range(x.getSize()):
        product *= math.exp(-x[i])
    return product

## We have to wrap f in a pysgpp.MultiFunction object.
func = pysgpp.multiFunc(f)

## Let's use a 3D-function.
d = 3


## @section combigrid_example_1 Example 1: Leja quadrature with linear growth of grid points
##
## Here comes the first and very simple example.
def example1():
    ## Let's increase the number of points by two for each level.
    growthFactor = 2

    ## Now create the operation object that handles the evaluation. The evaluation mode is quadrature,
    ## so it will approximate the integral of f over [0, 1]^d. It uses Leja points with 1 + 2*l
    ## points in level l. The level starts from zero, higher level means finer grid.
    ## Slower growth of the number of points per level means that the total number of points used can
    ## be controlled better.
    operation = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, func, growthFactor)

    ## Now, we compute the result. The parameter 2 means that grid at level-multi-indices with a
    ## 1-norm (i.e. sum of entries) less than or equal to 2 are used. In our 3D case, these are
    ## exactly the levels (0, 0, 0), (1, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 0), (0, 2, 0),
    ## (0, 1, 1), (0, 0, 1) and (0, 0, 2).
    result = operation.evaluate(2)

    ## Now compare the result to the analytical solution:
    print("Quadrature result: " + str(result) + ", analytical solution: " + str(math.pow(1.0 - 1.0/math.e, d)))

    ## We can also find out how many function evaluations have been used by accessing the storage
    ## which stores computed function values:
    print("Number of function evaluations: " + str(operation.getStorage().getNumEntries()))

## @section combigrid_example_2 Example 2: Polynomial interpolation on nested Clenshaw Curtis grids
##
## The next example uses interpolation.
def example2():
    ## This time, we use Clenshaw-Curtis points with exponentially growing number of points per level.
    ## This is helpful for CC points to make them nested. Nested means that the set of grid points at
    ## one level is a subset of the set of grid points at the next level. Nesting can drastically
    ## reduce the number of needed function evaluations. Using these grid points, we will do
    ## polynomial interpolation at a single point.
    operation = pysgpp.CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(d, func)

    ## Now create a point where to evaluate the interpolated function:
    evaluationPoint = pysgpp.DataVector([0.1572, 0.6627, 0.2378])

    ## We can now evaluate the interpolation at this point (using 3 as a bound for the 1-norm of the
    ## level multi-index):
    result = operation.evaluate(3, evaluationPoint)

    ## Now compare the result to the actual function value:
    print("Interpolation result: " + str(result) + ", function value: " + str(func(evaluationPoint)))

    ## Again, print the number of function evaluations:
    print("Function evaluations: " + str(operation.getStorage().getNumEntries()))

    ## Now, let's do another (more sophisticated) evaluation at a different point, so change the point
    ## and re-set the parameter. This method will automatically clear all intermediate values that
    ## have been computed internally up to now.
    evaluationPoint[0] = 0.4444
    print("Target function value: " + str(func(evaluationPoint)))
    operation.setParameters(evaluationPoint)

    ## The level manager provides more options for combigrid evaluation, so let's get it:
    levelManager = operation.getLevelManager()

    ## We can add regular levels like before:
    levelManager.addRegularLevels(3)

    ## The result can be fetched from the CombigridOperation:
    print("Regular result 1: " + str(operation.getResult()))
    print("Total function evaluations: " + str(operation.getStorage().getNumEntries()))

    ## We can also add more points in a regular structure, using at most 50 new function evaluations.
    ## All level-adding variants of levelManager also have a parallelized version. This version
    ## executes the calls to func in parallel with a specified number of threads, which is okay here
    ## since func supports parallel evaluations. Since func takes very little time to evaluate and the
    ## parallelization only concerns function evaluations and not the computations on the resulting
    ## function values, parallel evaluation is not actually useful in this case.
    ## We will use 4 threads for the function evaluations.
    levelManager.addRegularLevelsByNumPointsParallel(50, 4)
    print("Regular result 2: " + str(operation.getResult()))
    print("Total function evaluations: " + str(operation.getStorage().getNumEntries()))

    ## We can also use adaptive level generation. The adaption strategy depends on the subclass of
    ## LevelManager that is used. If you do not want to use the default LevelManager, you can specify
    ## your own LevelManager:
    operation.setLevelManager(pysgpp.AveragingLevelManager())
    levelManager = operation.getLevelManager()

    ## It was necessary to use setLevelManager(), because this links the LevelManager to the
    ## computation. Now, let's add at most 80 more function evaluations adaptively:
    levelManager.addLevelsAdaptive(60)
    print("Adaptive result: " + str(operation.getResult()))
    print("Total function evaluations: " + str(operation.getStorage().getNumEntries()))

## @section combigrid_example_3 Example 3: Evaluation at multiple points
##
## Now, we want to do interpolation at multiple evaluation points efficiently.
def example3():
    ## Use Leja points unlike example 2 and use CombigridMultiOperation for evaluation at multiple
    ## points.
    operation = pysgpp.CombigridMultiOperation.createLinearLejaPolynomialInterpolation(d, func)

    ## We slightly deviate from the C++ example here and pass the interpolation points via a DataMatrix.
    ## We will use 2 interpolation points.
    firstParam = [0.2, 0.6, 0.7]
    secondParam = [0.3, 0.9, 1.0]
    parameters = pysgpp.DataMatrix([firstParam, secondParam])

    ## Let's use the simple interface for this example and stop the time:
    stopwatch = pysgpp.Stopwatch()
    result = operation.evaluate(3, parameters)
    stopwatch.log()
    print("First result: " + str(result[0]) + ", function value: " + str(func(pysgpp.DataVector(firstParam))))
    print("Second result: " + str(result[1]) + ", function value: " + str(func(pysgpp.DataVector(secondParam))))

## @section combigrid_example_4 Example 4: Serialization and lookup tables
##
## This example shows how to store and retrieve computed function values.
## At first, we create a function that prints a string if it is called.
## This shows us when it is (not) called.
def loggingF(x):
    print("call function")
    return x[0]

def example4():
    ## After wrapping our new function into a pysgpp.MultiFunction, we create a FunctionLookupTable.
    ## This will cache the function values by their DataVector parameter and use cached values if available.
    ## Note, however, that even slightly differing DataVectors will lead to separate function evaluations.
    loggingFunc = pysgpp.multiFunc(loggingF)
    lookupTable = pysgpp.FunctionLookupTable(loggingFunc)
    operation = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, lookupTable.toMultiFunction())

    ## Do a normal computation...
    result = operation.evaluate(2)
    print("Result computed: " + str(result))

    ## The first (and most convenient) possibility to store the data is serializing the lookup table.
    ## The serialization is not compressed and will roughly use 60 Bytes per entry. If you have lots
    ## of data, you might consider compressing it.
    pysgpp.writeToFile("lookupTable.log", lookupTable.serialize())

    ## It is also possible to store which levels have been evaluated:
    pysgpp.writeToFile("levels.log", operation.getLevelManager().getSerializedLevelStructure())

    ## Restore the data into another lookup table. The function is still needed for new evaluations.
    restoredLookupTable = pysgpp.FunctionLookupTable(loggingFunc)
    restoredLookupTable.deserialize(pysgpp.readFromFile("lookupTable.log"))
    operation2 = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, restoredLookupTable.toMultiFunction())

    ## A new evaluation with the same levels does not require new function evaluations:
    operation2.getLevelManager().addLevelsFromSerializedStructure(pysgpp.readFromFile("levels.log"))
    result = operation2.getResult()
    print("Result computed (2nd time): " + str(result))

    ## Another less general way of storing the data is directly serializing the storage underlying the
    ## operation. This means that retrieval is faster, but it only works if the same grid is used
    ## again.
    ## For demonstration purposes, we use loggingFunc directly this time without a lookup table:
    pysgpp.writeToFile("storage.log", operation.getStorage().serialize())
    operation3 = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, pysgpp.multiFunc(loggingFunc))
    operation3.getStorage().deserialize(pysgpp.readFromFile("storage.log"))
    result = operation3.evaluate(2)
    print("Result computed (3rd time): " + str(result))

## @section combigrid_example_5 Example 5: Using different operations in each dimension
##
## This example shows how to apply different operators in different dimensions.
def example5():
    ## First, we want to configure which grid points to use in which dimension.
    ## We use Chebyshev points in the 0th dimension. To make them nested, we have to use at least \f$n
    ## = 3^l\f$ points at level \f$l\f$. This is why this method contains the prefix exp.
    ## CombiHierarchies provides some matching configurations for grid points. If you nevertheless
    ## need your own configuration or you want to know which growth strategy and ordering fit to which
    ## point distribution, look up the implementation details in CombiHierarchies, it is not
    ## difficult to implement your own configuration.
    grids = pysgpp.AbstractPointHierarchyVector()
    grids.push_back(pysgpp.CombiHierarchies.expChebyshev())

    ## Our next set of grid points are Leja points with linear growth (\f$n = 1 + 3l\f$).
    ## For the last dimension, we use equidistant points with boundary. These are suited for linear
    ## interpolation. To make them nested, again the slowest possible exponential growth is selected
    ## by the CombiHierarchies class.
    grids.push_back(pysgpp.CombiHierarchies.linearLeja(3))
    grids.push_back(pysgpp.CombiHierarchies.expUniform())

    ## The next thing we have to configure is the linear operation that is performed in those
    ## directions. We will use polynomial interpolation in the 0th dimension, quadrature in the 1st
    ## dimension and linear interpolation in the 2nd dimension.
    ## Roughly spoken, this means that a quadrature is performed on the 1D function that is the
    ## interpolated function with two fixed parameters. But since those operators "commute", the
    ## result is invariant under the order that the operations are applied in.
    ## The CombiEvaluators class also provides analogous methods and typedefs for the multi-evaluation
    ## case.
    evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
    evaluators.push_back(pysgpp.CombiEvaluators.polynomialInterpolation())
    evaluators.push_back(pysgpp.CombiEvaluators.quadrature())
    evaluators.push_back(pysgpp.CombiEvaluators.linearInterpolation())

    ## To create a CombigridOperation object with our own configuration, we have to provide a
    ## LevelManager as well:
    levelManager = pysgpp.WeightedRatioLevelManager()
    operation = pysgpp.CombigridOperation(grids, evaluators, levelManager, func)

    ## The two interpolations need a parameter \f$(x, z)\f$. If \f$\tilde{f}\f$ is the interpolated
    ## function, the operation approximates the result of \f$\int_0^1 \tilde{f}(x, y, z) \,dy\f$.
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
