from vectorFunctions import vectorObjFuncSGpp
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
from argparse import ArgumentParser
import sys
import os
import pickle
import time
matplotlib.use("TkAgg")


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


def execute(dim,
            numTimeSteps,
            gridType,
            degree,
            maxPoints,
            initialLevel,
            numRefine,
            stepType,
            numDraws,
            numTune,
            numChains,
            input_true,
            sigma_true,
            error_percent,
            numMeasurements):
    pyFunc = okushiri(dim, numTimeSteps=numTimeSteps)
    objFunc = vectorObjFuncSGpp(pyFunc)

    # set parameter ranges
    lb = pysgpp.DataVector(dim, 0.5)
    ub = pysgpp.DataVector(dim, 1.5)

    # create surrogate
    reSurf = pysgpp.SparseGridResponseSurfaceBsplineVector(
        objFunc, lb, ub, pysgpp.Grid.stringToGridType(gridType), degree)
    reSurf.surplusAdaptive(maxPoints, initialLevel, numRefine)
    print("created response surface with {} grid points".format(
        reSurf.getSize()))

    my_ODEop = ODEop(dim, reSurf, numTimeSteps)

    # create true measurement values
    trueValues = pysgpp.DataVector(numTimeSteps)
    objFunc.eval(input_true, trueValues)
    input_true_py = np.zeros(input_true.getSize())
    for d in range(input_true.getSize()):
        input_true_py[d] = input_true.get(d)
    # add errors to measurements.
    np.random.seed(0)
    measurements = np.zeros((numTimeSteps, numMeasurements))
    for m in range(numMeasurements):
        for t in range(numTimeSteps):
            measurements[t, m] = trueValues[t] + \
                (np.random.randn(1) * sigma_true)

    # The probabilistic model
    with pm.Model() as prob_model:

        # Priors for unknown model parameters
        all_params = pm.Uniform('parameter', lower=0.5, upper=1.5, shape=dim)

        # sigma = pm.HalfNormal('sigma', sd=0.1)
        sigma = pm.Lognormal('sigma', mu=-1, sd=1)

        # Forward model
        ode_sol = my_ODEop(all_params)
        forward = ode_sol

        # Likelihood
        Y_obs = []
        for m in range(numMeasurements):
            Y_obs.append(pm.Normal('Y_obs_%i' % m, mu=forward,
                                   sd=sigma, observed=measurements[:, m]))
        if stepType == 'Metropolis':
            step = pm.Metropolis()
        elif stepType == 'NUTS':
            step = pm.NUTS()
        start = time.time()
        trace = pm.sample(draws=numDraws, step=step, tune=numTune,
                          init='adapt_diag', chains=numChains)  # ,cores=4
        totalSampleTime = time.time()-start

        metaData = {
            'dim': dim,
            'numTimeSteps': numTimeSteps,
            'gridType': gridType,
            'degree': degree,
            'maxPoints': maxPoints,
            'initialLevel': initialLevel,
            'numRefine': numRefine,
            'numDraws': numDraws,
            'numTune': numTune,
            'numChains': numChains,
            'stepType': stepType,
            'numMeasurements': numMeasurements,
            'sigma_true': sigma_true,
            'input_true': input_true_py,
            'measurements': measurements,
            'error_percent': error_percent,
            'totalSampleTime': totalSampleTime
        }
        savePath = '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/traces'
        name = 'okushiri_{}D_{}_{}draws_{}tune_{}{}_{}grid_{}errp_{}msr'.format(
            dim, stepType, numDraws, numTune, gridType, degree, maxPoints, error_percent, numMeasurements)
        savePath = os.path.join(savePath, name)
        pm.save_trace(trace, directory=savePath, overwrite=True)
        with open(savePath+'/metaData.pkl', 'wb') as buff:
            pickle.dump(metaData, buff)
        gridStr = reSurf.serializeGrid()
        with open(savePath+'/grid.dat', 'w+') as fp:
            fp.write(gridStr)
        coeffs = reSurf.getCoefficients()
        coeffs.toFile(savePath + '/coeffs.dat')
        print('saved data to {}'.format(savePath))

        return trace, prob_model, measurements


############################ Main ############################
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    parser.add_argument('--demo', default='nothing', type=str,
                        help='this is how you add arguments')

    args = parser.parse_args()
    dim = 2
    numTimeSteps = 451
    gridType = 'nakbsplineboundary'
    degree = 3
    maxPoints = 500
    initialLevel = 1
    numRefine = 10

    stepType = 'NUTS'
    numDraws = 2000
    numTune = 1000
    numChains = 2
    numMeasurements = 10
    input_true = pysgpp.DataVector(dim, 0.77)
    # TODO: Check for the final domain we choose if the max is still 0.012!
    # Maximum true value is of size 0.012, p% error -> sigma = 0.01*p*0.012
    error_percent = 10
    sigma_true = 0.01*error_percent*0.012

    trace, prob_model, measurements = execute(dim, numTimeSteps, gridType, degree, maxPoints, initialLevel,
                                              numRefine, stepType, numDraws, numTune, numChains, input_true,
                                              sigma_true, error_percent, numMeasurements)

    print('Summary:')
    summary = pm.summary(trace)
    print(summary)

    print('\n')
    for d in range(dim):
        meanPosterior = summary['mean'][d]
        print('true parameter_{} - posterior mean = {}'.format(d,
                                                               input_true[d] - meanPosterior))
    print('true sigma - posterior mean        = {}'.format(
        sigma_true - summary['mean'][-1]))

    pm.traceplot(trace)
    plt.show()
