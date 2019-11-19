import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pysgpp
import warnings
from scipy.integrate import odeint
import time
import sys

# Bosch DC Motor
from dc_motor import dc_motor_analytical_I
from dc_motor import dc_motor_analytical_W
from dc_motor import dc_motor_ode_I
from dc_motor import dc_motor_ode_W

# ANUGA Okushiri Benchmark
sys.path.append('/home/rehmemk/git/anugasgpp/Okushiri')   # nopep8
from sgppOkushiri import okushiri, okushiri_g5, okushiri_g7, okushiri_g9, okushiri_input_wave  # nopep8
sys.path.append('/home/rehmemk/git/anugasgpp/Okushiri/investigation')   # nopep8
from sgppOkushiri_noResidual import okushiri_noResidual    # nopep8


def getFunction(model, dim=1, out=1, scalarModelParameter=16):
    if model == 'demo':
        return demo()

    # Bosch DC Motor
    elif model == 'dc_motor_analytical_I':
        return dc_motor_analytical_I(dim, numTimeSteps=out)
    elif model == 'dc_motor_analytical_W':
        return dc_motor_analytical_W(dim, numTimeSteps=out)
    elif model == 'dc_motor_ode_I':
        return dc_motor_ode_I(dim, numTimeSteps=out)
    elif model == 'dc_motor_ode_W':
        return dc_motor_ode_W(dim, numTimeSteps=out)
    elif model == 'lotkaVolterraHare':
        return lotkaVolterraHare(dim, numTimeSteps=out)

    # ANUGA Okushiri Benchmark
    elif model == 'okushiri':
        return okushiri(dim, numTimeSteps=out, gridResolution=scalarModelParameter)
    elif model == 'okushiri_g5':
        return okushiri_g5(dim, numTimeSteps=out, gridResolution=scalarModelParameter)
    elif model == 'okushiri_g7':
        return okushiri_g7(dim, numTimeSteps=out, gridResolution=scalarModelParameter)
    elif model == 'okushiri_g9':
        return okushiri_g9(dim, numTimeSteps=out, gridResolution=scalarModelParameter)
    elif model == 'okushiri_input_wave':
        return okushiri_input_wave(dim)
    elif model == 'maxOkushiri':
        return maxOkushiri(dim)
    elif model == 'okushiri_noResidual':
        return okushiri_noResidual(dim, numTimeSteps=out)

####################### auxiliary functions #######################


class vectorObjFuncSGpp(pysgpp.VectorFunction):
    # wraps the objective function for SGpp
    # NOTE: If we want to optimize we have to introduce an
    #       objFuncSGppSigned as in scalarFunction.py!

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

    def cleanUp(self):
        try:
            self.objFunc.cleanUp()
        except:
            warnings.warn('could not clean up')

    def getPrecalcData(self):
        return self.objFunc.getPrecalcData()


def randomPoint(lb, ub):
    dim = len(lb)
    point = np.random.rand(dim)
    for d in range(dim):
        point[d] = lb[d] + point[d] * (ub[d]-lb[d])
    return point


def createErrorDataSet(model, path, numMCPoints, dim=1, out=1, scalarModelParameter=3):
    # precalcultes evaluations of an objective function and stores them
    # so that they can be used for Monte carlo error calculations
    start = time.time()
    objFunc = getFunction(model, dim, out, scalarModelParameter)
    lb, ub = objFunc.getDomain()
    points = np.ndarray((numMCPoints, dim))
    evaluations = np.ndarray((out, numMCPoints))

    np.random.seed()
    for n in range(numMCPoints):
        sth = randomPoint(lb, ub)
        points[n, :] = sth
        this = objFunc.eval(points[n, :])
        evaluations[:, n] = this
        print("{} of {} points processed".format(n+1, numMCPoints))

    objFuncDataPath = os.path.join(path, objFunc.getName())
    if not os.path.exists(objFuncDataPath):
        os.makedirs(objFuncDataPath)
    pointsPath = os.path.join(
        objFuncDataPath, 'evaluationPoints.pkl')
    with open(pointsPath, 'wb') as fp:
        pickle.dump(points, fp)
        print('saved points to {}'.format(pointsPath))
    jacobianPath = os.path.join(
        objFuncDataPath, 'evaluations.pkl')
    with open(jacobianPath, 'wb') as fp:
        pickle.dump(evaluations, fp)
        print('saved evaluations to {}'.format(jacobianPath))
    print("precalculating {} values took {}s".format(
        numMCPoints, time.time()-start))


def createJacobianErrorDataSet(model, path, numMCPoints, dim=1, out=1, scalarModelParameter=3):
    # precalcultes evaluations of the jacobian of an objective function
    # and stores them so that they can be used for Monte carlo error calculations
    objFunc = getFunction(model, dim, out, scalarModelParameter)
    lb, ub = objFunc.getDomain()
    points = np.ndarray((numMCPoints, dim))
    jacobianEvaluations = np.ndarray((out, dim, numMCPoints))

    np.random.seed()
    for n in range(numMCPoints):
        sth = randomPoint(lb, ub)
        points[n, :] = sth
        this = objFunc.evalJacobian(points[n, :])
        jacobianEvaluations[:, :, n] = this
        print("{} of {} points processed".format(n+1, numMCPoints))

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


class lotkaVolterraHare():
    # The Lotka-Volterra model from the pymc3 example
    # https://docs.pymc.io/notebooks/ODE_parameter_estimation.html
    # Here we study how well the solution is approximated by our surrogate.
    # The model results in values for the hare and the lynx population.
    # This function only returns hare

    def __init__(self, dim, numTimeSteps=21):
        self.dim = dim
        self.out = numTimeSteps
        self.times = np.arange(0, numTimeSteps)

        # true values from the traceplot of the original example
        # extracted very inaccurate by eyesight, but here only used
        # to determine meaningful parameter ranges
        alpha_true = 0.55
        beta_true = 0.0265
        gamma_true = 0.775
        delta_true = 0.0225
        Xt0_true = 32.5
        Yt0_true = 5.8
        paramDeviation = 0.05
        lDev = 1 - paramDeviation
        rDev = 1 + paramDeviation

        self.pdfs = pysgpp.DistributionsVector()

        self.pdfs.push_back(pysgpp.DistributionUniform(0.25, 1.1))
        # lDev*alpha_true, rDev*alpha_true))
        if dim > 1:
            self.pdfs.push_back(pysgpp.DistributionUniform(0.01, 0.06))
            # lDev*beta_true, rDev*beta_true))
        if dim > 2:
            self.pdfs.push_back(pysgpp.DistributionUniform(0.45, 1.35))
            # lDev*gamma_true, rDev*gamma_true))
        if dim > 3:
            self.pdfs.push_back(pysgpp.DistributionUniform(0.01, 0.06))
            # lDev*delta_true, rDev*delta_true))
        if dim > 4:
            self.pdfs.push_back(pysgpp.DistributionUniform(9.0, 47.5))
            # lDev*Xt0_true, rDev*Xt0_true))
        if dim > 5:
            self.pdfs.push_back(pysgpp.DistributionUniform(4.0, 11.0))
            # lDev*Yt0_true, rDev*Yt0_true))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub

    def getName(self):
        return "lotkaVolterraHare{}D{}T".format(self.getDim(), self.getOut())

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.out

    def eval(self, x):
        alpha = x[0]
        beta = 0.0265
        gamma = 0.775
        delta = 0.0225
        Xt0 = 32.5
        Yt0 = 5.8
        if self. dim > 1:
            beta = x[1]
        if self.dim > 2:
            gamma = x[2]
        if self.dim > 3:
            delta = x[3]
        if self.dim > 4:
            Xt0 = x[4]
        if self.dim > 5:
            Yt0 = x[5]
        parameters = [alpha, beta, gamma, delta, Xt0, Yt0]

        def r(y, t, p):
            X, Y = y
            dX_dt = alpha*X - beta*X*Y
            dY_dt = -gamma*Y + delta*X*Y
            return dX_dt, dY_dt

        # The rtol and atol parameters determine the quality of the ODE solution
        # By default they are ca 1e-8.
        # With this the B-spline surrogate error stagnates at ca 1000 points at
        # around 1e-6. If we want a better surrogate with more grid points,
        # smaller rtol and atol have to be used.
        # A quick check showed 50% increase in runtime when using 1e-10 instead
        # of default 1e-8
        values = odeint(r, [Xt0, Yt0], self.times,
                        (parameters,))  # , rtol=1e-10, atol=1e-10)
        hare = values[:, 0]
        lynx = values[:, 1]
        return hare

    def evalJacobian(self, x):
        print('Jacobian unknown')
        return 0

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
