import matplotlib.pyplot as plt
import pymc3 as pm
import numpy as np
import matplotlib
from collections import OrderedDict
matplotlib.use("TkAgg")


def histogram(trace, val, true_value, meanPosterior, sdPosterior, numBins=50, lineplot=1, barplot=0):
    numChains = trace.nchains

    if len(val) % numChains == 0:
        chainLength = int(len(val)/numChains)
    else:
        print('Warning: The number of chains ({}) and values ({}) don\'t match'.format(
            numChains, len(val)))

    chainColors = ['C9', 'C4', 'C6', 'C2']
    for c in range(numChains):
        chain_val = val[c*chainLength:(c+1)*chainLength]
        hist, bin_edges = np.histogram(chain_val, bins=numBins, density=True)
        if lineplot == 1:
            left, right = bin_edges[:-1], bin_edges[1:]
            X = np.array([left, right]).T.flatten()
            Y = np.array([hist, hist]).T.flatten()
            plt.plot(
                X, Y, '-', color=chainColors[c % len(chainColors)], label='chain {}'.format(c))
    # add the maximum likelihood normal distribution
    mu = meanPosterior
    sigma = sdPosterior
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bin_edges - mu))**2))
    plt.plot(bin_edges, y, color='C1', linewidth=2,
             label='posterior distribution')

    # add 5 and 95 percentiles calculated over ALL chains
    fifth, ninetyfifth = np.percentile(val, [5, 95])
    plt.plot(np.ones(2)*fifth, [0, np.max(hist)],
             'C7', label='5%, 95% percentiles')
    plt.plot(np.ones(2)*ninetyfifth, [0, np.max(hist)], 'C7')

    # plot true value and posterior mean
    plt.plot(np.ones(2)*true_value, [0, np.max(hist)], 'k', label='true value')
    plt.plot(np.ones(2)*meanPosterior,
             [0, np.max(hist)], 'C3', label='posterior mean')
    # plot posterior standard deviation
    plt.plot([meanPosterior-0.5*sdPosterior, meanPosterior+0.5*sdPosterior],
             [np.max(hist)/4, np.max(hist)/4], 'C0', label='posterior sigma')
    plt.plot([meanPosterior-0.5*sdPosterior, meanPosterior-0.5*sdPosterior],
             [np.max(hist)/4 - np.max(hist)/50, np.max(hist)/4 + np.max(hist)/50], 'C0')
    plt.plot([meanPosterior+0.5*sdPosterior, meanPosterior+0.5*sdPosterior],
             [np.max(hist)/4 - np.max(hist)/50, np.max(hist)/4 + np.max(hist)/50], 'C0')

    if barplot == 1:
        width = 0.7 * (bin_edges[1] - bin_edges[0])
        center = (bin_edges[:-1] + bin_edges[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
    plt.gca().get_yaxis().set_visible(False)
    handles, labels = plt.gca().get_legend_handles_labels()
    return handles, labels


def ppc(trace,
        prob_model,
        ppcMeasurementIndex,
        forward,
        sigma,
        measurements,
        input_true,
        pyFunc,
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

        # This plots ALL measurements at once
        # so numMeasurements*numTimeSteps many measurement points
        plt.plot(range(len(measurements)), measurements[:, 0],
                 '+', color='b', label='Observed')
        plt.plot(range(len(measurements)), mean_ppc,
                 color='g', marker='*', label='mean of ppc')
        plt.plot(range(len(measurements)), CriL_ppc, '--', color='b',
                 lw=2, label='credible intervals')
        plt.plot(range(len(measurements)), CriU_ppc, '--', color='b', lw=2)

        trueValues = pyFunc.eval([input_true[0], input_true[1]])
        plt.plot(range(len(trueValues)), trueValues,
                 color='r', label='true vallues')

        # every measurement series gets the same label
        # remove multiple labels in this step
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
