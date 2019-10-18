import sys
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
from collections import OrderedDict
from sgppOkushiri import okushiri
from okushiri_Bayes import ODEop
matplotlib.use("TkAgg")


sys.path.append('/home/rehmemk/git/SGpp/MR_Python/Bayes')  # nopep8
from plotRoutinesBayes import ppc  # nopep8
from plotRoutinesBayes import histogram  # nopep8


def oneHistogram(trace, metaData, qoi, numBins=50, lineplot=1, barplot=0):
    # qoi: 0,1,2,3 for the paramters, -1 for sigma
    if qoi >= 0:
        val = trace.get_values('parameter')[:, qoi]
        true_value = metaData['input_true'][qoi]
    elif qoi == -1:
        val = trace.get_values('sigma')
        true_value = metaData['sigma_true']
    else:
        print('Warning: qoi {} does not exist'.format(qoi))
    meanPosterior = summary['mean'][qoi]
    sdPosterior = summary['sd'][qoi]

    handles, labels = histogram(trace, val, true_value, meanPosterior, sdPosterior,
                                numBins, lineplot, barplot)

    return handles, labels


def allHistograms(trace, metaData, numBins=50, lineplot=1, barplot=0):
    dim = metaData['dim']
    usedVariables = list(range(dim))
    usedVariables.append(-1)
    fig = plt.figure()
    for i, qoi in enumerate(usedVariables):
        plt.subplot(np.ceil(len(usedVariables)/2.), 2, i+1)
        if qoi >= 0:
            plt.title('parameter {}'.format(qoi))
        elif qoi == -1:
            plt.title('sigma')
        handles, labels = oneHistogram(
            trace, metaData, qoi, numBins, lineplot, barplot)
    fig.legend(handles, labels, loc='lower center', ncol=3)
    plt.tight_layout()


def loadTrace(dim, gridType, degree, maxPoints, numDraws, numTune, stepType, error_percent, numMeasurements):
    path = '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/traces'
    name = 'okushiri_{}D_{}_{}draws_{}tune_{}{}_{}grid_{}errp_{}msr'.format(
        dim, stepType, numDraws, numTune, gridType, degree, maxPoints, error_percent, numMeasurements)
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
error_percent = 10
numMeasurements = 10

trace, metaData, prob_model, reSurf, forward, sigma = loadTrace(dim, gridType,
                                                                degree, maxPoints, numDraws, numTune, stepType, error_percent, numMeasurements)

measurements = metaData['measurements']
numTimeSteps = metaData['numTimeSteps']
input_true = metaData['input_true']
sigma_true = metaData['sigma_true']
summary = pm.summary(trace)
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


do_ppc = 1
if do_ppc == 1:
    ppcMeasurementIndex = 0
    pyFunc = okushiri(dim, numTimeSteps=numTimeSteps)
    ppc(trace, prob_model, ppcMeasurementIndex, forward,
        sigma, measurements, input_true, pyFunc)

do_histogram = 0
if do_histogram == 1:
    qoi = 0
    plt.figure()
    histogram(trace, metaData, qoi)

do_allHistograms = 1
if do_allHistograms == 1:
    allHistograms(trace, metaData)

do_pymc3_traceplot = 0
if do_pymc3_traceplot == 1:
    pm.traceplot(trace)

plt.show()
