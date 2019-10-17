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

def histogram(trace, metaData, qoi, numBins=50, lineplot=1, barplot=0, legend = 1):
    # qoi: 0,1,2,3 for the paramters, -1 for sigma
    if qoi >= 0:
        val = trace.get_values('parameter')[:,qoi]
        true_value = metaData['input_true'][qoi]
    elif qoi == -1:
        val = trace.get_values('sigma')
        true_value = metaData['sigma_true']
    else:
        print('Warning: qoi {} does not exist'.format(qoi))
    meanPosterior = summary['mean'][qoi]
    sdPosterior = summary['sd'][qoi]

    numChains = trace.nchains
    if len(val) % numChains == 0:
        chainLength = int(len(val)/numChains)
    else:
        print('Warning: The number of chains and values don\'t match')

    chainColors = ['C9', 'C4', 'C6', 'C2']
    for c in range(numChains):
        chain_val = val[c*chainLength:(c+1)*chainLength]
        hist, bin_edges = np.histogram(chain_val, bins=numBins, density=True)
        if lineplot == 1:
            left, right = bin_edges[:-1], bin_edges[1:]
            X = np.array([left, right]).T.flatten()
            Y = np.array([hist, hist]).T.flatten()
            plt.plot(X, Y, '--', color=chainColors[c%len(chainColors)], label='chain {}'.format(c))
    # add the maximum likelihood normal distribution
    mu = meanPosterior
    sigma = sdPosterior
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bin_edges - mu))**2))
    plt.plot(bin_edges, y, color='C1', linewidth=2, label='distribution')
    # plot true value and posterior mean
    plt.plot(np.ones(2)*true_value, [0, np.max(hist)], 'k', label='true value')
    plt.plot(np.ones(2)*meanPosterior, [0, np.max(hist)], 'C3', label='posterior mean')
    # plot posterior standard deviation
    plt.plot([meanPosterior-0.5*sdPosterior, meanPosterior+0.5*sdPosterior], [np.max(hist)/4, np.max(hist)/4], 'C0', label='posterior sigma')
    plt.plot([meanPosterior-0.5*sdPosterior, meanPosterior-0.5*sdPosterior], [np.max(hist)/4 - np.max(hist)/50, np.max(hist)/4 + np.max(hist)/50], 'C0')
    plt.plot([meanPosterior+0.5*sdPosterior, meanPosterior+0.5*sdPosterior], [np.max(hist)/4 - np.max(hist)/50, np.max(hist)/4 + np.max(hist)/50], 'C0')

    if barplot == 1:
        width = 0.7 * (bin_edges[1] - bin_edges[0])
        center = (bin_edges[:-1] + bin_edges[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
    plt.gca().get_yaxis().set_visible(False)  
    if legend == 1:
        plt.legend()
    else:
        handles, labels = plt.gca().get_legend_handles_labels()
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
        handles, labels = histogram(trace, metaData, qoi, numBins, lineplot, barplot, legend = 0)
    fig.legend(handles, labels, loc='lower center', ncol=3)
    plt.tight_layout()


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
    ppc(prob_model, ppcMeasurementIndex, forward,
        sigma, measurements, numTimeSteps, input_true)

do_histogram = 0
if do_histogram == 1:
    qoi = 0
    plt.figure()
    histogram(trace, metaData, qoi)

do_allHistograms = 1
if do_allHistograms == 1:
    allHistograms(trace, metaData)

do_pymc3_traceplot = 1
if do_pymc3_traceplot == 1:
    pm.traceplot(trace)

plt.show()
