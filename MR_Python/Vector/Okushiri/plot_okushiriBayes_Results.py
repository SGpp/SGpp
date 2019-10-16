from vectorFunctions import vectorObjFuncSGpp
import matplotlib.pyplot as plt
import pymc3 as pm
import numpy as np
import theano
import theano.tensor as T
from scipy.integrate import odeint
import pysgpp
import pickle
import os
import matplotlib
from sgppOkushiri import okushiri
from okushiri_Bayes import ODEop
matplotlib.use("TkAgg")


def ppc(prob_model,
        ppcMeasurementIndex,
        forward,
        sigma,
        measurements,
        numTimeSteps,
        input_true,
        numSamples=500):
    # draw posterior predictives
    # that is sample the posterior parameter distributions 'numSamples' many times
    # then calculate the okushiri benchmark for each sample resulting in 'numSamples'
    # many timelines
    # TODO: Look into Y_obs parameter. currently uses the measuremnt given by
    # ppcMeasurementIndex, what about the rest?
    with prob_model:
        pm.Normal('Y_obs', mu=forward, sd=sigma,
                  observed=measurements[:, ppcMeasurementIndex])

        ppc_samples = pm.sample_posterior_predictive(
            trace, numSamples, model=prob_model)['Y_obs']
        # we plot the mean and two percentiles of the ppc data and the measured data
        # this demonstrates that our posterior guesses for the parameters
        # lead to okushiri outputs fitting the data well.
        mean_ppc = ppc_samples.mean(axis=0)
        # the two percentiles cover 95% of the samples
        CriL_ppc = np.percentile(ppc_samples, q=2.5, axis=0)
        CriU_ppc = np.percentile(ppc_samples, q=97.5, axis=0)

        plt.figure()
        plt.plot(range(len(measurements)), measurements,
                 '+', color='b', label='Observed')
        plt.plot(range(len(measurements)), mean_ppc,
                 color='g', marker='*', label='mean of ppc')
        plt.plot(range(len(measurements)), CriL_ppc, '--', color='b',
                 lw=2, label='credible intervals')
        plt.plot(range(len(measurements)), CriU_ppc, '--', color='b', lw=2)

        pyFunc = okushiri(dim, numTimeSteps=numTimeSteps)
        trueValues = pyFunc.eval([input_true[0], input_true[1]])
        plt.plot(range(len(trueValues)), trueValues,
                 color='r', label='true vallues')
        plt.legend()


def loadTrace(dim, gridType, degree, maxPoints, numDraws, numTune, stepType):
    path = '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/traces'
    name = 'okushiri_{}D_{}_{}draws_{}tune_{}{}_{}'.format(
        dim, stepType, numDraws, numTune, gridType, degree, maxPoints)
    path = os.path.join(path, name)

    with open(path+'/metaData.pkl', 'rb') as buff:
        metaData = pickle.load(buff)
    dim = metaData['dim']
    numTimeSteps = metaData['numTimeSteps']
    degree = metaData['degree']
    with open(path+'/grid.dat', 'r') as fp:
        gridStr = fp.read()
    coeffs = pysgpp.DataMatrix_fromFile(path + '/coeffs.dat')

    pyFunc = okushiri(dim, numTimeSteps=numTimeSteps)
    objFunc = vectorObjFuncSGpp(pyFunc)
    lb = pysgpp.DataVector(dim, 0.5)
    ub = pysgpp.DataVector(dim, 1.5)
    reSurf = pysgpp.SparseGridResponseSurfaceBsplineVector(
        objFunc, lb, ub, gridStr, degree, coeffs)
    my_ODEop = ODEop(dim, reSurf, numTimeSteps)
    with pm.Model() as prob_model:

        # Priors for unknown model parameters
        all_params = pm.Uniform('parameter', lower=0.5,
                                upper=1.5, shape=dim)
        # sigma = pm.HalfNormal('sigma', sd=0.1)
        sigma = pm.Lognormal('sigma', mu=-1, sd=1)

        ode_sol = my_ODEop(all_params)
        forward = ode_sol

        trace = pm.load_trace(path)
        return trace, metaData, prob_model, reSurf, forward, sigma


dim = 2
gridType = 'nakbsplineboundary'
degree = 3
maxPoints = 500
numDraws = 2000
numTune = 1000
stepType = 'NUTS'

trace, metaData, prob_model, reSurf, forward, sigma = loadTrace(dim, gridType,
                                                                degree, maxPoints, numDraws, numTune, stepType)

summary = pm.summary(trace)
print(summary)
pm.traceplot(trace)

measurements = metaData['measurements']
numTimeSteps = metaData['numTimeSteps']
input_true = metaData['input_true']
do_ppc = 1
if do_ppc == 1:
    ppcMeasurementIndex = 0
    ppc(prob_model, ppcMeasurementIndex, forward,
        sigma, measurements, numTimeSteps, input_true)

plt.show()
