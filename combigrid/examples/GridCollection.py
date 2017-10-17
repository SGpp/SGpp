import pysgpp
import FunctionClasses as fctClass
import matplotlib.pyplot as plt
import numpy as np
import random
import sympy as sp


class CombiCombigrid1d_hermite:
    def __init__(self, function, grad_function):
        self.func = pysgpp.multiFunc(function)
        self.grad_func = pysgpp.multiFunc(grad_function)
        self.d = 1
        self.operation_psi = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)
        self.operation_zeta = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
            self.d, self.grad_func)

    def evaluate(self, level, x):
        return self.operation_psi.evaluate(level,
                                           x) + 1 * self.operation_zeta.evaluate(
            level, x)

    def getLevelManager(self):
        return self.operation_psi.getLevelManager()


class CombiCombigrid2dHermite_without_mixed:
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.d = dim

        self.grad = []
        for i in range(self.d):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.operation_psi = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)

        self.operation_zeta = []
        for i in range(self.d):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                    self.d, i, self.grad[i]))

    def evaluate(self, level, x):
        sum = 0
        sum += self.operation_psi.evaluate(level, x)

        for op in self.operation_zeta:
            sum += op.evaluate(level, x)

        return sum

    def getLevelManager(self):
        return self.operation_psi.getLevelManager()


class CombiCombigridHermite:
    def __init__(self, func_container, dim):
        self.d = dim
        self.func_container = func_container

        self.func = pysgpp.multiFunc(self.func_container.getFunction())
        self.grad = []

        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)

        # create all combinations for mixed gradients
        self.mixed_directions = list(fctClass.powerset([i for i in range(self.d)]))
        self.mixed_grad = []
        self.operationzeta_psi = []

        for i in range(1, len(self.mixed_directions)):
            mixed_grad_temp = self.func_container.getGradient(self.mixed_directions[i])
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))
        for i in range(len(self.mixed_grad)):
            self.operationzeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaInterpolation(self.d,
                                                                                    self.mixed_directions[
                                                                                        i + 1],
                                                                                    self.mixed_grad[
                                                                                        i]))

    def evaluate(self, level, x):
        sum = 0

        sum += self.operation_psi.evaluate(level, x)
        # for op in self.operation_zeta_psi:
        #    sum += op.evaluate(level, x)


        for op in self.operationzeta_psi:
            sum += op.evaluate(level, x)
        return sum

    def getLevelManager(self):
        return self.operation_psi.getLevelManager()


class BaseLinearFullgrid():
    def __init__(self, dim, func_container, full=False):
        self.dim = dim

        self.func = function_eval_undisplaced(func_container.getFunction())
        self.level = -1
        self.full = full

    def linearGrid(self, level):
        grid = pysgpp.Grid.createLinearBoundaryGrid(self.dim)
        self.gridStorage = grid.getStorage()

        if (self.full == False):
            grid.getGenerator().regular(level)
        else:
            grid.getGenerator().full(level)
        alpha = pysgpp.DataVector(self.gridStorage.getSize())
        alpha.setAll(0.0)
        for i in range(self.gridStorage.getSize()):
            gp = self.gridStorage.getPoint(i)
            coordinates = pysgpp.DataVector(self.dim)
            gp.getStandardCoordinates(coordinates)
            alpha[i] = self.func.evalUndisplaced(coordinates)

        pysgpp.createOperationHierarchisation(grid).doHierarchisation(alpha)

        self.alpha = alpha
        return grid

    def evaluate(self, level, x):
        if self.level != level:
            self.level = level
            self.grid = self.linearGrid(self.level)
            self.opEvalLinear = pysgpp.createOperationEvalNaive(self.grid)

        return self.opEvalLinear.eval(self.alpha, x)

    def gridOperationEval(self, opEval, alpha):
        def eval(x):
            return opEval.eval(alpha, x)

        return eval

    def getLevelManager(self):
        return self.LevelManager(self.gridStorage)

    class LevelManager:
        def __init__(self, gridstorage):
            self.gridStorage = gridstorage

        def numGridPoints(self):
            return self.gridStorage.getSize()

    class function_eval_undisplaced:
        def __init__(self, function):
            self.function = function

        def evalUndisplaced(self, x):
            return self.function(x)


class HierachGridBSpline:
    # create a bsplineGrid with calculated hierarchized coefficients

    def __init__(self, dim, degree, func_collection, full=False):
        self.dim = dim
        self.degree = degree
        self.func = function_eval_undisplaced(func_collection.getFunction())
        self.level = -1
        self.full = full

    def createBsplineGrid(self, level):

        grid = pysgpp.Grid.createBsplineBoundaryGrid(self.dim, self.degree)
        self.gridStorage = grid.getStorage()
        print "dimensionality:         {}".format(self.gridStorage.getDimension())

        if (self.full == False):
            grid.getGenerator().regular(level)
        else:
            grid.getGenerator().full(level)

        functionValues = pysgpp.DataVector(self.gridStorage.getSize())
        functionValues.setAll(0.0)
        for i in xrange(self.gridStorage.getSize()):
            gp = self.gridStorage.getPoint(i)
            coordinates = pysgpp.DataVector(self.dim)
            gp.getStandardCoordinates(coordinates)
            # print( coordinates[0],coordinates[1])
            functionValues[i] = self.func.evalUndisplaced(coordinates)

        print ("Hierarchizing...\n")
        coeffs = pysgpp.DataVector(len(functionValues))
        hierSLE = pysgpp.OptHierarchisationSLE(grid)
        sleSolver = pysgpp.OptAutoSLESolver()

        if not sleSolver.solve(hierSLE, functionValues, coeffs):
            print "Solving failed, exiting."

        print(self.level)

        return grid, coeffs

    def gridOperationEval(self, opEval, alpha):
        def eval(x):
            return opEval.eval(alpha, x)

        return eval

    def evaluate(self, level, x):
        if self.level != level:
            self.level = level
            self.grid, self.coeffs = self.createBsplineGrid(self.level)
            self.opEvalBspline = pysgpp.createOperationEvalNaive(self.grid)

        return self.opEvalBspline.eval(self.coeffs, x)

    # gets created on function call, no "refresh"
    def getLevelManager(self):
        return self.LevelManager(self.gridStorage, self.dim)

    class LevelManager:
        def __init__(self, gridstorage, dim):
            self.gridStorage = gridstorage
            self.dim = dim

        def numGridPoints(self):
            return self.gridStorage.getSize()

        def getAllGridPoints(self):
            list = []
            for i in xrange(self.gridStorage.getSize()):
                gp = self.gridStorage.getPoint(i)
                coordinates = pysgpp.DataVector(self.dim)
                gp.getStandardCoordinates(coordinates)
                list.append(coordinates)

            return list


class function_eval_undisplaced:
    def __init__(self, function):
        self.function = function

    def evalUndisplaced(self, x):
        return self.function(x)


class CombiCombigriddeBaarHarding:
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_psi = []
        self.operation_zeta = []
        for i in range(dim):
            self.operation_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(self.d,
                                                                                         i,
                                                                                         self.func))

        for i in range(dim):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
                    self.d, i, self.grad[i]))

    def evaluate(self, level, x):
        sum = 0

        for operation in self.operation_psi:
            sum += operation.evaluate(level, x)
        for operation in self.operation_zeta:
            sum += operation.evaluate(level, x)
        sum -= self.operation_linear.evaluate(level, x) * (self.d - 1)
        return sum

    def getLevelManager(self):
        return self.operation_linear.getLevelManager()


class CombiCombigriddeBaarHardingBSpline:

    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_psi = []
        self.operation_zeta = []
        for i in range(dim):
            self.operation_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(self.d,
                                                                                         i,
                                                                                         self.func))

        for i in range(dim):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
                    self.d, i, self.grad[i]))

    def evaluate(self, level, x):
        sum = 0

        for operation in self.operation_psi:
            sum += operation.evaluate(level, x)
        for operation in self.operation_zeta:
            sum += operation.evaluate(level, x)
        sum -= self.operation_linear.evaluate(level, x) * (self.d - 1)
        return sum

    def getLevelManager(self):
        return self.operation_linear.getLevelManager()


################FullGrids########################################################
class FullGrid:
    def eval_full(self, level, x, fullgrid_operation):
        scalars = []
        for i in range(len(x)):
            scalars.append(pysgpp.FloatScalarVector(x[i]))

        scalars = pysgpp.FloatScalarVectorVector(scalars)
        fullgrid_operation.setParameters(scalars)

        return fullgrid_operation.eval([level for i in range(self.d)]).getValue()

    class LevelManager:
        def __init__(self, operation, level, d):
            self.operation = operation
            self.level = level
            self.d = d

        # returns points with fullgrid
        def numGridPoints(self):
            return self.operation.numPoints([self.level for i in range(self.d)])

        def getAllGridPoints(self):
            print self.level
            return self.operation.getGridPoints([self.level for i in range(self.d)])


class CombiFullGridHermite(FullGrid):
    def __init__(self, func_collection, dim):
        self.d = dim
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        self.func_collection = func_collection

        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)
        self.operation_psi = self.operation_psi.getFullGridEvaluator()

        # create all combinations for the mixed gradients
        self.mixed_directions = list(fctClass.powerset([i for i in range(self.d)]))
        self.mixed_grad = []
        self.operationzeta_psi = []

        for i in range(1, len(self.mixed_directions)):
            mixed_grad_temp = self.func_collection.getGradient(self.mixed_directions[i])
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))

        for i in range(len(self.mixed_grad)):
            self.operationzeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaInterpolation(self.d,
                                                                                    self.mixed_directions[
                                                                                        i + 1],
                                                                                    self.mixed_grad[
                                                                                        i]))
        # lazy full conversion to full grid operations

        for i in range(len(self.operationzeta_psi)):
            self.operationzeta_psi[i] = self.operationzeta_psi[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level

        sum += self.eval_full(level, x, self.operation_psi)
        # for op in self.operation_zeta_psi:
        #    sum += op.evaluate(level, x)


        for op in self.operationzeta_psi:
            sum += self.eval_full(level, x, op)
        return sum

    def getLevelManager(self):
        return self.LevelManager(self.operation_psi, self.level, self.d)


class CombiFullGridHermite_withoutmixed(FullGrid):
    def __init__(self, func_collection, dim):
        self.d = dim
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        self.func_collection = func_collection

        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)
        self.operation_psi = self.operation_psi.getFullGridEvaluator()

        # create all combinations for the mixed gradients

        self.mixed_grad = []
        self.operationzeta_psi = []

        for i in range(self.d):
            mixed_grad_temp = self.func_collection.getGradient([i])
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))

        for i in range(self.d):
            self.operationzeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaInterpolation(self.d,
                                                                                    [i],
                                                                                    self.mixed_grad[
                                                                                        i]))
        # lazy full conversion to full grid operations

        for i in range(len(self.operationzeta_psi)):
            self.operationzeta_psi[i] = self.operationzeta_psi[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level

        sum += self.eval_full(level, x, self.operation_psi)
        # for op in self.operation_zeta_psi:
        #    sum += op.evaluate(level, x)


        for op in self.operationzeta_psi:
            sum += self.eval_full(level, x, op)
        return sum

    def getLevelManager(self):
        return self.LevelManager(self.operation_psi, self.level, self.d)


class CombiFullGridLinear(FullGrid):
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_linear = self.operation_linear.getFullGridEvaluator()

        self.operation_psi = []
        self.operation_zeta = []
        for i in range(dim):
            self.operation_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(self.d,
                                                                                         i,
                                                                                         self.func))
        # lazy redefine of fullgrid
        for i in range(dim):
            self.operation_psi[i] = self.operation_psi[i].getFullGridEvaluator()

        for i in range(dim):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
                    self.d, i, self.grad[i]))

        # lazy redefine of fullgrid
        for i in range(dim):
            self.operation_zeta[i] = self.operation_zeta[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level

        for operation in self.operation_psi:
            sum += self.eval_full(level, x, operation)
        for operation in self.operation_zeta:
            sum += self.eval_full(level, x, operation)
        sum -= self.eval_full(level, x, self.operation_linear) * (self.d - 1)
        return sum

    def getLevelManager(self):
        return self.LevelManager(self.operation_linear, self.level, self.d)


class LinearFullgrid(FullGrid):
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_linear = self.operation_linear.getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level
        return self.eval_full(level, x, self.operation_linear)

    def getLevelManager(self):
        return self.LevelManager(self.operation_linear, self.level, self.d)
