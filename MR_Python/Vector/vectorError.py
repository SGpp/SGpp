import numpy as np
import pysgpp
import sys


def l2Vector(reSurfVector, objFunc, numMCPoints):
    dim = objFunc.getDim()
    out = objFunc.getOut()
    lb = reSurfVector.getLowerBounds()
    ub = reSurfVector.getUpperBounds()

    outwiseErrors = np.zeros(out)
    trueEval = pysgpp.DataVector(dim)
    for _ in range(numMCPoints):
        point = np.random.rand(dim)
        point = [lb[d]+(ub[d]-lb[d])*point[d] for d in range(dim)]
        objFunc.eval(point, trueEval)
        reSurfEval = reSurfVector.eval(pysgpp.DataVector(point))
        for d in range(dim):
            outwiseErrors[d] += (trueEval[d]-reSurfEval[d])**2
    outwiseErrors = np.sqrt(outwiseErrors/numMCPoints)
    return outwiseErrors


def nrmseVector(reSurfVector, objFunc, numMCPoints):
    dim = objFunc.getDim()
    out = objFunc.getOut()
    lb = reSurfVector.getLowerBounds()
    ub = reSurfVector.getUpperBounds()

    outwiseErrors = np.zeros(out)
    func_max = np.ones(out)*(-sys.float_info.max)
    func_min = np.ones(out)*sys.float_info.max
    trueEval = pysgpp.DataVector(dim)
    for _ in range(numMCPoints):
        point = np.random.rand(dim)
        point = [lb[d]+(ub[d]-lb[d])*point[d] for d in range(dim)]
        objFunc.eval(point, trueEval)
        reSurfEval = reSurfVector.eval(pysgpp.DataVector(point))
        for d in range(dim):
            outwiseErrors[d] += (trueEval[d]-reSurfEval[d])**2
            if trueEval[d] > func_max[d]:
                func_max[d] = trueEval[d]
            if trueEval[d] < func_min[d]:
                func_min[d] = trueEval[d]
    outwiseErrors = np.sqrt(outwiseErrors/numMCPoints)
    for d in range(d):
        outwiseErrors[d] /= (func_max[d]-func_min[d])
    return outwiseErrors, func_min, func_max
