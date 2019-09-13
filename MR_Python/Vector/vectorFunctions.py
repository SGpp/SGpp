import os
import pickle
import numpy as np
import pysgpp
import matplotlib.pyplot as plt

# MALTE: Er nimmt als schwierige Funktion immer sin(sin(x))
# Extrem schwierig: sin(1/x), oszilliert enorm am linken Rand


# Functions are evaluated in a point given as DataVector.
def getFunction(model, dim=1, out=1, scalarModelParameter=3):
    if model == 'demo':
        return demo()

####################### auxiliary functions #######################
# wraps the objective function for SGpp
# NOTE: If we want to optimize we have to introduce an objFuncSGppSigned as in scalarFunction.py!


class vectorObjFuncSGpp(pysgpp.VectorFunction):

    # input dimension dim
    # output dimension out
    def __init__(self, objFunc):
        self.dim = objFunc.getDim()
        self.out = objFunc.getOut()
        self.objFunc = objFunc
        super().__init__(self.dim, self.out)

    def eval(self, x, value):
        result = self.objFunc.eval(x)
        for t in range(self.out):
            value.set(t, result[t])

    def evalJacobian(self, x):
        jacobian = self.objFunc.evalJacobian(x)
        return jacobian

    def getName(self):
        return self.objFunc.getName()

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.out

    def getLowerBounds(self):
        lb, _ = self.objFunc.getDomain()
        return lb

    def getUpperBounds(self):
        _, ub = self.objFunc.getDomain()
        return ub

    def getDistributions(self):
        return self.objFunc.getDistributions()

    def getMean(self):
        return self.objFunc.getMean()

    def getVar(self):
        return self.objFunc.getVar()


def randomPoint(lb, ub):
    dim = len(lb)
    point = np.random.rand(dim)
    for d in range(dim):
        point[d] = lb[d] + point[d] * (ub[d]-lb[d])
    return point


def createJacobianErrorDataSet(model, path, numMCPoints, dim=1, out=1, scalarModelParameter=3):
    objFunc = getFunction(model, dim, out, scalarModelParameter)
    lb, ub = objFunc.getDomain()
    # print(lb.toString())
    # print(ub.toString())
    points = np.ndarray((numMCPoints, dim))
    jacobianEvaluations = np.ndarray((out, dim, numMCPoints))

    np.random.seed()
    for n in range(numMCPoints):
        sth = randomPoint(lb, ub)
        points[n, :] = sth
        this = objFunc.evalJacobian(points[n, :])
        jacobianEvaluations[:, :, n] = this

        # print(points[n, :])
        # print(jacobianEvaluations[:, :, n])

    objFuncDataPath = os.path.join(path, objFunc.getName())
    if not os.path.exists(objFuncDataPath):
        os.makedirs(objFuncDataPath)
    pointsPath = os.path.join(
        objFuncDataPath, 'evaluationPoints.pkl')
    with open(pointsPath, 'wb') as fp:
        pickle.dump(points, fp)
        print('saved points to {}'.format(pointsPath))
    jacobianPath = os.path.join(
        objFuncDataPath, 'jacobianEvaluations.pkl')
    with open(jacobianPath, 'wb') as fp:
        pickle.dump(jacobianEvaluations, fp)
        print('saved jacobian evaluations to {}'.format(jacobianPath))

#########################################################
#####################   functions   #####################
#########################################################


class demo():
    # this demonstrates what a fully equipped function looks like
    # IT also serves as a test function for developing/testing

    def __init__(self):
        self.dim = 2
        self.out = 3
        self.pdfs = pysgpp.DistributionsVector()

        self.pdfs.push_back(pysgpp.DistributionUniform(-2, 2))
        self.pdfs.push_back(pysgpp.DistributionNormal(0.0, 0.1))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub

    def getName(self):
        return "demo{}D{}O".format(self.getDim(), self.getOut())

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.out

    def eval(self, x):
        result = np.zeros(self.getOut())
        result[0] = x[0] + 2*x[1]
        result[1] = (x[0]-x[1])**3
        result[2] = np.exp(x[0]+2*x[1])
        return result

    def evalJacobian(self, x):
        jacobian = np.zeros((self.getOut(), self.getDim()))
        jacobian[0, 0] = 1
        jacobian[1, 0] = 3*(x[0]-x[1])**2
        jacobian[2, 0] = np.exp(x[0]+2*x[1])

        jacobian[0, 1] = 2
        jacobian[1, 1] = -3*(x[0]-x[1])**2
        jacobian[2, 1] = 2*np.exp(x[0]+2*x[1])
        return jacobian

    def getDistributions(self):
        return self.pdfs

    def getMeans(self):
        print("Means unknown")
        means = np.zeros(self.getOut())
        return means

    def getVars(self):
        print("Vars unknown")
        variances = np.zeros(self.getOut())
        return variances

    def getIntegrals(self):
        print("Integrals unknown")
        integrals = np.zeros(self.getOut())
        return integrals
