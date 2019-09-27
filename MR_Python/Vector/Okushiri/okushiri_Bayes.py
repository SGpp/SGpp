import ipdb
import matplotlib
import pysgpp
from scipy.integrate import odeint
import theano.tensor as T
import theano
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
from sgppOkushiri import okushiri
import sys
matplotlib.use("TkAgg")

# TODO: This should be included from vectorFunctions.py!


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
        self.objFunc.cleanUp()


##########
dim = 2
numTimeSteps = 451
gridType = 'nakbsplineboundary'
degree = 3
maxPoints = 500
initialLevel = 1
numRefine = 10
numDraws = 2000
numTune = 1000
numChains = 2
##########

pyFunc = okushiri(dim, numTimeSteps=numTimeSteps)
objFunc = vectorObjFuncSGpp(pyFunc)

# set parameter ranges
lb = pysgpp.DataVector(dim, 0.0)
ub = pysgpp.DataVector(dim, 1.0)

# create surrogate
reSurf = pysgpp.SparseGridResponseSurfaceBsplineVector(
    objFunc, lb, ub, pysgpp.Grid.stringToGridType(gridType), degree)
reSurf.surplusAdaptive(maxPoints, initialLevel, numRefine)
print("created I response surface with {} grid points".format(
    reSurf.getSize()))

# Testing one eval
# somePoint = pysgpp.DataVector(dim, 0.33)
# someEval = reSurf.eval(somePoint)
# print("Here is something for you:")
# print(someEval)

# measure error of surrogates
# numErrPoints = 1000
# componentwiseError = pysgpp.DataVector(len(times))
# averageL2Error = reSurf.averageL2Error(
#     I_motor_instance, IComponentwiseError, numErrPoints)
# print("average I error: {}    (min {} max {})".format(
#     IAverageL2Error[0], IAverageL2Error[1], IAverageL2Error[2]))


class ODEGradop(theano.Op):
    # auxiliary class for gradient evaluations
    def __init__(self, dim, surrogate, numTimeSteps):
        self.dim = dim
        self.surrogate = surrogate
        self.numTimeSteps = numTimeSteps

    def make_node(self, x, g):
        x = theano.tensor.as_tensor_variable(x)
        g = theano.tensor.as_tensor_variable(g)
        node = theano.Apply(self, [x, g], [g.type()])
        return node

    def perform(self, node, inputs_storage, output_storage):
        x = inputs_storage[0]
        g = inputs_storage[1]
        out = output_storage[0]

        x_sgpp = pysgpp.DataVector(x)
        jacobian = pysgpp.DataMatrix(self.numTimeSteps, self.dim)
        _ = self.surrogate.evalJacobian(x_sgpp, jacobian)

        # # TODO: I changed this from a 3D matrix to a 2D matrix.
        # # Is it still correctly shaped and filled?
        sensitivities = np.zeros((self.numTimeSteps, 1, self.dim))
        for t in range(self.numTimeSteps):
            for d in range(dim):
                sensitivities[t, 0, d] = jacobian.get(t, d)
        sensitivities = sensitivities.reshape(
            (1*self.numTimeSteps, len(x)))
        out[0] = sensitivities.T.dot(g)


class ODEop(theano.Op):
    # class to evaluate the surrogates and their gradients (via auxiliary class)

    def __init__(self, dim, surrogate, numTimeSteps):
        self.dim = dim
        self.surrogate = surrogate
        self.numTimeSteps = numTimeSteps

    def make_node(self, x):
        x = theano.tensor.as_tensor_variable(x)
        return theano.Apply(self, [x], [x.type()])

    def perform(self, node, inputs_storage, output_storage):
        x = inputs_storage[0]
        out = output_storage[0]

        x_sgpp = pysgpp.DataVector(x)
        evaluation = self.surrogate.eval(x_sgpp)
        State = np.ones((evaluation.getSize(), ))
        for t in range(evaluation.getSize()):
            State[t] = evaluation[t]
        out[0] = State

    def grad(self, inputs, output_grads):
        x = inputs[0]
        g = output_grads[0]
        grad_op = ODEGradop(self.dim, self.surrogate, self.numTimeSteps)
        grad_op_apply = grad_op(x, g)
        return [grad_op_apply]


my_ODEop = ODEop(dim, reSurf, numTimeSteps)

# create true measurement values
trueValues = pysgpp.DataVector(numTimeSteps)
input_true = pysgpp.DataVector(dim, 0.33)
objFunc.eval(input_true, trueValues)
# add errors to measurements. Maximum true value is of size 0.012, 1% error
sigma_true = 0.00012  # =0.01*0.012
np.random.seed(0)
Y = np.zeros((numTimeSteps, ))
for t in range(numTimeSteps):
    Y[t] = trueValues[t] + (np.random.randn(1) * sigma_true)

# The probabilistic model
with pm.Model() as prob_model:

    # Priors for unknown model parameters
    param1 = pm.Uniform('parameter 1', lower=0.0, upper=1.0)
    param2 = pm.Uniform('parameter 2', lower=0.0, upper=1.0)
    all_params = pm.math.stack([param1, param2], axis=0)

    # sigma = pm.HalfNormal('sigma', sd=0.1)
    sigma = pm.Lognormal('sigma', mu=-1, sd=1)

    # Forward model
    ode_sol = my_ODEop(all_params)
    forward = ode_sol

    # Likelihood
    Y_obs = pm.Normal('Y_obs', mu=forward, sd=sigma, observed=Y)
    # Y_obs = pm.Lognormal('Y_obs', mu=pm.math.log(
    #     forward), sd=sigma, observed=Y)
    step = pm.Metropolis()
    trace = pm.sample(draws=numDraws, step=step, tune=numTune,
                      init='adapt_diag', chains=numChains)  # ,cores=4

with prob_model:
    pm.traceplot(trace)
plt.show()
