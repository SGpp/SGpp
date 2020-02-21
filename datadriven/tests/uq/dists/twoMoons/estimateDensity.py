# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
import sys
import os
import json
import csv
import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
from scipy.special import binom

from pysgpp import KernelType_EPANECHNIKOV, KernelType_GAUSSIAN
from pysgpp.extensions.datadriven.uq.dists import SGDEdist, KDEDist
from pysgpp.extensions.datadriven.uq.plot import plotSG2d
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature
from pysgpp.extensions.datadriven.tools import writeGrid, writeAlphaARFF
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d
from pysgpp.extensions.datadriven.uq.dists.Lognormal import Lognormal
from pysgpp.extensions.datadriven.uq.dists.Beta import Beta
from pysgpp.extensions.datadriven.uq.dists.J import J
from pysgpp.pysgpp_swig import BandwidthOptimizationType_MAXIMUMLIKELIHOOD, \
    createOperationMakePositive, DataVector, \
    BandwidthOptimizationType_SILVERMANSRULE
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d

from positivity.plot import plotResults, plotCostsPerIteration
from argparse import ArgumentParser
from positivity.helper_functions import strToCandidateSearchAlgorithm, \
    strToInterpolationAlgorithm, strTobandwidthOptimizationType
from scipy.stats.stats import kstest
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform
from pysgpp.extensions.datadriven.uq.dists.MultivariateNormal import MultivariateNormal
from pysgpp.extensions.datadriven.uq.dists.Dist import Dist
from pysgpp.extensions.datadriven.uq.dists.NatafDist import NatafDist


def splitset(samples, splitPercentage=0.8):
    numSamples = samples.shape[0]
    numTrainSamples = int(np.ceil(numSamples * splitPercentage))
    trainSamplesixs = np.random.choice(numSamples, numTrainSamples, replace=False)

    trainSamplesMask = np.zeros(numSamples, dtype="bool")
    trainSamplesMask[trainSamplesixs] = True
    testSamplesMask = np.invert(trainSamplesMask)

    assert numTrainSamples == samples[trainSamplesMask].shape[0]
    
    return samples[trainSamplesMask], samples[testSamplesMask]


def estimateSGDEDensity(functionName,
                        trainSamples,
                        testSamples=None,
                        bounds=None,
                        iteration=0,
                        plot=False,
                        out=True,
                        label="sgde_zero",
                        candidates="intersections", interpolation="setToZero"):
    print("train: %i x %i (mean=%g, var=%g)" % (trainSamples.shape[0], trainSamples.shape[1], np.mean(trainSamples), np.var(trainSamples)))
    if testSamples is not None:
        print("test : %i x %i (mean=%g, var=%g)" % (testSamples.shape[0], testSamples.shape[1], np.mean(testSamples), np.var(testSamples)))

    candidateSearchAlgorithm = strToCandidateSearchAlgorithm(candidates)
    interpolationAlgorithm = strToInterpolationAlgorithm(interpolation)

    results = {}
    crossEntropies = {}
    config = {"grid_level": 1,
              "grid_type": "linear",
              "grid_maxDegree": 1,
              "refinement_numSteps": 0,
              "refinement_numPoints": 3,
              "solver_threshold": 1e-10,
              "solver_verbose": False,
              "regularization_type": "Laplace",
              "crossValidation_enable": True,
              "crossValidation_kfold": 5,
              "crossValidation_silent": True,
              "sgde_makePositive": False}

    pathResults = os.path.join("data", label)
    key = 1
    bestCV = float("Inf")
    bestDist = None

    # stats
    stats = {'config': {'functionName': functionName,
                        'numDims': 2,
                        'adaptive': True,
                        'refnums': 0,
                        'consistentGrid': True,
                        'candidateSearchAlgorithm': candidates,
                        'interpolationAlgorithm': interpolation,
                        'maxNumGridPoints': 0,
                        'iteration': iteration},
             'trainSamples': trainSamples,
             'testSamples': testSamples}

    for level in range(2, 7):
        print("-" * 60)
        print("l=%i" % level)
        for refinementSteps in range(0, 5):
            config["grid_level"] = level
            config["refinement_numSteps"] = refinementSteps
            sgdeDist = SGDEdist.byLearnerSGDEConfig(trainSamples, config=config,
                                                    bounds=bounds)
            # -----------------------------------------------------------
            grid, alpha = sgdeDist.grid, sgdeDist.alpha
            cvSgde = sgdeDist.crossEntropy(testSamples)

            maxLevel = grid.getStorage().getMaxLevel()
            numDims = grid.getStorage().getDimension()

            print("  " + "-" * 30)
            print("  #ref = %i: gs=%i -> CV test = %g" % (refinementSteps,
                                                          sgdeDist.grid.getSize(),
                                                          cvSgde))
            # -----------------------------------------------------------
            # make it positive
            positiveGrid = grid.clone()
            positiveAlpha_vec = DataVector(alpha)
            opPositive = createOperationMakePositive(candidateSearchAlgorithm,
                                                     interpolationAlgorithm,
                                                     True, False)
            opPositive.makePositive(positiveGrid, positiveAlpha_vec, True)

            # scale to unit integrand
            positiveAlpha = positiveAlpha_vec.array()
            positiveSgdeDist = SGDEdist(positiveGrid, positiveAlpha, trainSamples,
                                        bounds=bounds)
            # -----------------------------------------------------------
            cvPositiveSgde = positiveSgdeDist.crossEntropy(testSamples)

            if plot and numDims == 2:
                fig = plt.figure()
                plotSG2d(grid, alpha, show_negative=True, show_grid_points=True)
                plt.title("pos: N=%i: vol=%g, log=%g" % (positiveGrid.getSize(),
                                                         doQuadrature(positiveGrid, positiveAlpha),
                                                         cvPositiveSgde))
                plt.tight_layout()
                if out:
                    plt.savefig(os.path.join(pathResults, "%s_density_pos_i%i_l%i_r%i.jpg" % (label, iteration, level, refinementSteps)))
                    plt.savefig(os.path.join(pathResults, "%s_density_pos_i%i_l%i_r%i.pdf" % (label, iteration, level, refinementSteps)))
                else:
                    plt.close(fig)


            # -----------------------------------------------------------
            print("  positive: gs=%i -> CV test = %g" % (positiveGrid.getSize(), cvPositiveSgde))
            # -----------------------------------------------------------
            # select the best density available based on the given criterion
            results[key] = {'config': config,
                            'dist': positiveSgdeDist}
            crossEntropies[key] = cvPositiveSgde
            key += 1
            candidateSearch = opPositive.getCandidateSetAlgorithm()
     
            if cvPositiveSgde < bestCV:
                bestCV = cvPositiveSgde
                bestDist = positiveSgdeDist 
                numComparisons = candidateSearch.costsComputingCandidates()

                # update the stats -> just for the current best one
                # write the stats of the current best results to the stats dict
                C = np.ndarray(numDims - 1, dtype="int")
                M = np.sum([1 for i in range(len(alpha)) if alpha[i] < 0])
                for d in range(2, numDims + 1):
                    C[d - 2] = binom(M, d)

                stats['config']['refnums'] = refinementSteps
                stats['config']['adaptive'] = refinementSteps > 0
                stats['negSGDE_json'] = sgdeDist.toJson()
                stats['posSGDE_json'] = positiveSgdeDist.toJson()
                stats['level'] = level
                stats['maxLevel'] = maxLevel
                stats['fullGridSize'] = (2 ** maxLevel - 1) ** numDims
                stats['sparseGridSize'] = grid.getSize()
                stats['discretizedGridSize'] = positiveGrid.getSize()
                stats['crossEntropyTrainZeroSGDE'] = sgdeDist.crossEntropy(trainSamples)
                stats['crossEntropyTrainDiscretizedSGDE'] = positiveSgdeDist.crossEntropy(trainSamples)
                stats['crossEntropyTestZeroSGDE'] = cvSgde
                stats['crossEntropyTestDiscretizedSGDE'] = cvPositiveSgde
                stats['numCandidates'] = int(candidateSearch.numCandidates())
                stats['numCandidatesPerLevel'] = np.array(candidateSearch.numCandidatesPerLevel().array(), dtype="int")
                stats['numCandidatesPerIteration'] = np.array(candidateSearch.numCandidatesPerIteration().array(), dtype="int")
                stats['costsCandidateSearch'] = candidateSearch.costsComputingCandidates()
                stats['costsCandidateSearchBinomial'] = int(C.sum())
                stats['costsCandidateSearchPerIteration'] = np.array(candidateSearch.costsComputingCandidatesPerIteration().array(), dtype="int")
                stats['costsCandidateSearchPerIterationBinomial'] = C

                if plot and numDims == 2:
                    fig = plt.figure()
                    plotSG2d(positiveGrid, positiveAlpha, show_negative=True, show_grid_points=False,
                             colorbarLabel=r"$f_{\mathcal{I}^\text{SG} \cup \mathcal{I}^\text{ext}}$")
                    plt.title(r"positive: $N=%i/%i$; \# comparisons$=%i$" % (positiveGrid.getSize(),
                                                                             (2 ** maxLevel - 1) ** numDims,
                                                                             numComparisons))
                    plt.xlabel(r"$\xi_1$")
                    plt.ylabel(r"$\xi_2$")
#                     plt.title(r"N=%i $\rightarrow$ %i: log=%g $\rightarrow$ %g" % (sgdeDist.grid.getSize(),
#                                                                                    positiveSgdeDist.grid.getSize(),
#                                                                                    cvSgde,
#                                                                                    cvPositiveSgde))
                    plt.tight_layout()
                    plt.savefig(os.path.join(pathResults, "%s_pos_i%i_l%i_r%i.jpg" % (label, iteration, level, refinementSteps)))
                    plt.savefig(os.path.join(pathResults, "%s_pos_i%i_l%i_r%i.pdf" % (label, iteration, level, refinementSteps)))
                    if out:
                        plt.close(fig)

                    fig, ax, _ = plotSG3d(positiveGrid, positiveAlpha)
                    ax.set_zlabel(r"$f_{\mathcal{I}^{\text{SG}} \cup \mathcal{I}^\text{ext}}(\xi_1, \xi_2)$", fontsize=20)
                    ax.set_xlabel(r"$\xi_1$", fontsize=20)
                    ax.set_ylabel(r"$\xi_2$", fontsize=20)

                    plt.tight_layout()
                    plt.savefig(os.path.join(pathResults, "%s_pos_i%i_l%i_r%i_3d.jpg" % (label, iteration, level, refinementSteps)))
                    plt.savefig(os.path.join(pathResults, "%s_pos_i%i_l%i_r%i_3d.pdf" % (label, iteration, level, refinementSteps)))
                    if out:
                        plt.close(fig)

            if plot and numDims == 2 and not out:
                plt.show()


    if out:
        # save stats
        filename = os.path.join("data", label, "stats_d%i_a%i_r%i_i%i_%s_%s.pkl" % (numDims, 1, refinementSteps, iteration, candidates, interpolation))
        fd = open(filename, "w")
        pkl.dump(stats, fd)
        fd.close()
        print("stats saved to -> '%s'" % filename)

        # dictionary that stores the information on the estimated densities
        myjson = {"Grid": {"dimNames": ["phi", "log(K_A)"],
                           "matrixEntries": ["phi", "log(K_A)"]},
                  "Set": {"path": "",
                          "grids": [],
                          "alphas": [],
                          "paramValues": [],
                          "paramName": "grid_size"}}
    
        for key, result in list(results.items()):
            config = result['config']
            dist = result['dist']
            # serialize grid and coefficients
            out = "sgde.i%i.k%i.N%i" % (iteration, key, dist.grid.getSize())
            out_grid = os.path.join(pathResults, "%s.grid" % out)
            out_alpha = os.path.join(pathResults, "%s.alpha.arff" % out)
            writeGrid(out_grid, dist.grid)
            writeAlphaARFF(out_alpha, dist.alpha)
    
            # collect information for json
            myjson["Set"]["grids"].append(os.path.abspath(out_grid))
            myjson["Set"]["alphas"].append(os.path.abspath(out_alpha))
            myjson["Set"]["paramValues"].append(crossEntropies[key])
            # -----------------------------------------------------------
            # serialize the config
            out_config = os.path.join(pathResults, "sgde.i%i.k%i.config" % (iteration, key))
            fd = open(out_config, "w")
            json.dump(config, fd, ensure_ascii=True, indent=True)
            fd.close()
    
            crossEntropies[key] = (crossEntropies[key], out_grid, out_alpha, out_config)
    
        # sort the results in myjson according to the cross entropy
        ixs = np.argsort(myjson["Set"]["paramValues"])
        myjson["Set"]["grids"] = [myjson["Set"]["grids"][ix] for ix in ixs]
        myjson["Set"]["alphas"] = [myjson["Set"]["alphas"][ix] for ix in ixs]
        myjson["Set"]["paramValues"] = [myjson["Set"]["paramValues"][ix] for ix in ixs]
    
        # serialize myjson
        out_config = os.path.join(pathResults, "sgde_visualization.i%i.config" % iteration)
        fd = open(out_config, "w")
        json.dump(myjson, fd, ensure_ascii=True, indent=True)
        fd.close()
    
        # serialize cross entropies
        out_crossEntropies = os.path.join(pathResults, "sgde_cross_entropies.i%i.csv" % iteration)
        fd = open(out_crossEntropies, 'wb')
        file_writer = csv.writer(fd)
        file_writer.writerow(["crossEntropy", "grid", "alpha", "sgdeConfig"])
        for out in list(crossEntropies.values()):
            file_writer.writerow(out)
        fd.close()
    
        # serialize samples
        np.savetxt(os.path.join(pathResults, "sgde_train_samples.i%i.csv" % iteration), trainSamples)
        np.savetxt(os.path.join(pathResults, "sgde_test_samples.i%i.csv" % iteration), testSamples)

        # serialize best configuration to json
        out_bestDist = os.path.join(pathResults, "sgde_best_config.i%i.json" % iteration)
        text = bestDist.toJson()
        fd = open(out_bestDist, "w")
        fd.write(text)
        fd.close()

    return bestDist, stats


def estimateKDEDensity(functionName,
                       trainSamples,
                       testSamples=None,
                       iteration=0,
                       plot=False,
                       out=True,
                       label="kde_gaussian",
                       bandwidthOptimizationTypeStr="rot"):
    print("train: %i x %i (mean=%g, var=%g)" % (trainSamples.shape[0], trainSamples.shape[1], np.mean(trainSamples), np.var(trainSamples)))
    if testSamples is not None:
        print("test : %i x %i (mean=%g, var=%g)" % (testSamples.shape[0], testSamples.shape[1], np.mean(testSamples), np.var(testSamples)))

    if "gaussian" in label:
        kernelType = KernelType_GAUSSIAN
    elif "epanechnikov" in label:
        kernelType = KernelType_EPANECHNIKOV
    else:
        raise AttributeError("label is unknown")

    bandwidthOptimizationType = strTobandwidthOptimizationType(bandwidthOptimizationTypeStr)
    kdeDist = KDEDist(trainSamples,
                      kernelType=kernelType,
                      bandwidthOptimizationType=bandwidthOptimizationType)
    # -----------------------------------------------------------
    cvKDE = kdeDist.crossEntropy(testSamples)

    if plot and kdeDist.getDim() == 2:
        fig = plt.figure()
        plotDensity2d(kdeDist)
        plt.title("log=%g" % cvKDE)
        if out:
            plt.tight_layout()
            plt.savefig(os.path.join(pathResults, "kde_dist.%s.i%i.jpg" % (functionName, iteration)))
            plt.savefig(os.path.join(pathResults, "kde_dist.%s.i%i.pdf" % (functionName, iteration)))
            if out:
                plt.close(fig)
        else:
            plt.show()

    print("CV test = %g" % cvKDE)

    # -----------------------------------------------------------
    if out:        
        pathResults = os.path.join("data", label)
    
        # serialize cross entropies
        out_crossEntropies = os.path.join(pathResults, "kde_cross_entropies.%s.i%i.csv" % (functionName, iteration))
        fd = open(out_crossEntropies, 'wb')
        file_writer = csv.writer(fd)
        file_writer.writerow(["crossEntropy"])
        file_writer.writerow([cvKDE])
        fd.close()
    
        # serialize samples
        np.savetxt(os.path.join(pathResults, "kde_train_samples.%s.i%i.csv" % (functionName, iteration)), trainSamples)
        np.savetxt(os.path.join(pathResults, "kde_test_samples.%s.i%i.csv" % (functionName, iteration)), testSamples)

        if plot:
            # plot density
            fig = plt.figure()
            plotDensity2d(kdeDist)
            plt.title("%s -> CV = %g" % (kdeDist.getBandwidths(), cvKDE))
            plt.savefig(os.path.join(pathResults, "kde_pdf.%s.i%i.jpg" % (functionName, iteration)))
            plt.close(fig)

        # serialize best configuration to json
        out_bestDist = os.path.join(pathResults, "kde_best_config.%s.i%i.json" % (functionName, iteration))
        text = kdeDist.toJson()
        fd = open(out_bestDist, "w")
        fd.write(text)
        fd.close()

    # stats
    stats = {'config': {'functionName': functionName,
                        'numDims': 2,
                        'label': label,
                        'bandwidth_optimization': BandwidthOptimizationType_MAXIMUMLIKELIHOOD, 
                        'kernelType': kernelType,
                        'iteration': iteration},
             'trainSamples': trainSamples,
             'testSamples': testSamples,
             'crossEntropyTrainKDE': kdeDist.crossEntropy(trainSamples),
             'crossEntropyTestKDE': cvKDE,
             'KDEDist_json': kdeDist.toJson()}

    return kdeDist, stats


def estimateNatafDensity(functionName,
                         natafType,
                         testSamples=None,
                         iteration=0,
                         bounds=None,
                         plot=False,
                         out=True,
                         label="nataf"):
    if "samples" in natafType:
        trainSamples = natafType["samples"]
        print("train: %i x %i (mean=%g, var=%g)" % (trainSamples.shape[0], trainSamples.shape[1], np.mean(trainSamples), np.var(trainSamples)))
    if testSamples is not None:
        print("test : %i x %i (mean=%g, var=%g)" % (testSamples.shape[0], testSamples.shape[1], np.mean(testSamples), np.var(testSamples)))

    # -----------------------------------------------------------
    if natafType["name"] == "samples":
        natafDist = NatafDist.by_samples(natafType["samples"].T)
    elif natafType["name"] == "gamma":
        natafDist = NatafDist.gamma_marginals(natafType["alpha"],
                                              natafType["beta"],
                                              natafType["cov"],
                                              bounds)
    elif natafType["name"] == "beta":
        natafDist = NatafDist.beta_marginals(natafType["lwr"],
                                             natafType["upr"],
                                             natafType["alpha"],
                                             natafType["beta"],
                                             natafType["cov"],
                                             bounds)
    elif natafType["name"] == "normal":
        natafDist = NatafDist.normal_marginals(natafType["mean"],
                                               natafType["stddev"],
                                               natafType["cov"],
                                               bounds)
    else:
        raise AttributeError("nataf type '%s' is not supported" % natafType["name"])

    cvNataf = natafDist.crossEntropy(testSamples)

    if plot and natafDist.getDim() == 2:
        fig = plt.figure()
        plotDensity2d(natafDist)
        plt.title("log=%g" % cvNataf)
        if out:
            plt.tight_layout()
            plt.savefig(os.path.join(pathResults, "nataf_dist.%s.i%i.jpg" % (functionName, iteration)))
            plt.savefig(os.path.join(pathResults, "nataf_dist.%s.i%i.pdf" % (functionName, iteration)))
            if out:
                plt.close(fig)
        else:
            plt.show()

    print("CV test = %g" % cvNataf)

    # -----------------------------------------------------------
    if out:
        pathResults = os.path.join("data", label)

        # serialize cross entropies
        out_crossEntropies = os.path.join(pathResults, "nataf_cross_entropies.%s.i%i.csv" % (functionName, iteration))
        fd = open(out_crossEntropies, 'wb')
        file_writer = csv.writer(fd)
        file_writer.writerow(["crossEntropy"])
        file_writer.writerow([cvNataf])
        fd.close()

        # serialize samples
        if "samples" in natafType:
            np.savetxt(os.path.join(pathResults, "nataf_train_samples.%s.i%i.csv" % (functionName, iteration)), trainSamples)
        np.savetxt(os.path.join(pathResults, "nataf_test_samples.%s.i%i.csv" % (functionName, iteration)), testSamples)

        if plot:
            # plot density
            fig = plt.figure()
            plotDensity2d(natafDist)
            plt.title("CV = %g" % (cvNataf,))
            plt.savefig(os.path.join(pathResults, "nataf_pdf.%s.i%i.jpg" % (functionName, iteration)))
            plt.close(fig)

        # serialize best configuration to json
        out_bestDist = os.path.join(pathResults, "nataf_best_config.%s.i%i.json" % (functionName, iteration))
        text = natafDist.toJson()
        fd = open(out_bestDist, "w")
        fd.write(text)
        fd.close()

    # stats
    stats = {'config': {'functionName': functionName,
                        'numDims': 2,
                        'label': label,
                        'cov': natafDist.cov(),
                        'iteration': iteration},
             'testSamples': testSamples,
             'crossEntropyTestNataf': cvNataf,
             'NatafDist_json': natafDist.toJson()}

    return natafDist, stats



def load_data_set(data_set, numSamples, numDims=2):
    np.random.seed(1234567)
    natafType = {}

    if "mult" in data_set:
        corr = 0.005
        var = 0.01
        diag = np.diag(np.ones(numDims) * var)
        offdiag = (np.ones((numDims, numDims)) - np.diag(np.ones(numDims))) * corr
        covMatrix = diag + offdiag

        natafType["cov"] = covMatrix
        if "normal" in data_set:
            mean = 0.5
            stddev = np.sqrt(covMatrix[0, 0])
            U = MultivariateNormal(np.ones(numDims) * mean, covMatrix, 0, 1)
            samples = U.rvs(numSamples)
            bounds = U.getBounds()
            natafType["name"] = "normal"
            natafType["mean"] = mean
            natafType["stddev"] = stddev

        elif "gamma" in data_set:
            alpha = 2
            beta = 3
            bounds = np.array([[0, 20]] * numDims)
            U = NatafDist.gamma_marginals(alpha, beta, covMatrix,
                                          bounds=bounds)
            rvs_samples = U.rvs(numSamples)
            # remove all samples which are outside the domain
            samples = np.array([])
            for sample in rvs_samples:
                isValid = True
                for idim in range(numDims):
                    isValid &= np.all(sample[idim] < bounds[0][1])
                if isValid:
                    samples = np.append(samples, sample)
            samples = samples.reshape((samples.size // numDims), numDims)

            natafType["name"] = "gamma"
            natafType["alpha"] = alpha
            natafType["beta"] = beta
        elif "beta" in data_set:
            alpha = 5.
            beta = 10.
            bounds = np.array([[0, 1]] * numDims)
            U = NatafDist.beta_marginals(0, 1, alpha, beta, covMatrix,
                                         bounds=bounds)
            rvs_samples = U.rvs(numSamples)
            # remove all samples which are outside the domain
            samples = np.array([])
            for sample in rvs_samples:
                isValid = True
                for idim in range(numDims):
                    isValid &= np.all(sample[idim] < bounds[0][1])
                if isValid:
                    samples = np.append(samples, sample)
            samples = samples.reshape((samples.size // numDims), numDims)

            natafType["name"] = "beta"
            natafType["alpha"] = alpha
            natafType["beta"] = beta
            natafType["lwr"] = 0
            natafType["upr"] = 1
    else:
        if "moons" in data_set:
            samples = np.loadtxt("data/twomoons.csv")
            bounds = np.array([[0.0, 1.0], [0.0, 1.0]])
        elif "friedman" in data_set:
            samples = np.loadtxt("data/friedman2_4d_50000.csv")
            bounds = None
        else:
            raise AttributeError()
        if samples.shape[0] > numSamples:
            ixs = np.random.randint(0, samples.shape[0], numSamples)
            samples = samples[ixs, :]
        natafType = {"name": "samples",
                     "samples": samples}

    return samples, bounds, natafType


def run_densityEstimation(functionName,
                          method,
                          kfold=20,
                          numDims=2,
                          numSamples=1000,
                          candidates="join",
                          bandwidthOptimizationType=BandwidthOptimizationType_SILVERMANSRULE,
                          out=True,
                          plot=False,
                          tikz=False):
    if method == "sgde_zero":
        interpolation = "zero"
    else:  # interpolation == "boundaries":
        interpolation = "boundaries"

    samples, bounds, natafType = load_data_set(functionName, numSamples, numDims)

    # do kfold cross validation
    crossEntropyValidation = np.zeros((kfold, 2))
    learnSamples, validationSamples = splitset(samples, splitPercentage=0.7)

    stats = {}
    for i in range(kfold):
        print("=" * 100)
        print("run (%s)= %i/%i" % (method, i + 1, kfold))
        print("=" * 100)
        print("valid: %i x %i (mean=%g, var=%g)" % (validationSamples.shape[0], validationSamples.shape[1], np.mean(validationSamples), np.var(validationSamples)))

        np.random.seed(i * 123456 + i % 2)
        trainSamples, testSamples = splitset(learnSamples, splitPercentage=1. - 1. / kfold)

        if "sgde" in method:
            dist, stats[i] = estimateSGDEDensity(functionName,
                                                 trainSamples,
                                                 testSamples,
                                                 bounds=bounds,
                                                 iteration=i,
                                                 plot=plot,
                                                 label=method,
                                                 out=out,
                                                 candidates=candidates,
                                                 interpolation=interpolation)
        elif "kde" in method:
            dist, stats[i] = estimateKDEDensity(functionName,
                                                trainSamples,
                                                testSamples,
                                                iteration=i,
                                                plot=plot,
                                                label=method,
                                                out=out,
                                                bandwidthOptimizationTypeStr=bandwidthOptimizationType)
        elif "nataf" in method:
            # estimate nataf density
            dist, stats[i] = estimateNatafDensity(functionName,
                                                  natafType,
                                                  testSamples,
                                                  iteration=i,
                                                  bounds=bounds,
                                                  plot=plot,
                                                  label=method,
                                                  out=out)
        else:
            raise AttributeError("unknown config '%s'" % method)

        # evaluate the distribution according to the validation set
        crossEntropyValidation[i, 0] = i
        crossEntropyValidation[i, 1] = dist.crossEntropy(validationSamples)
        stats[i]["crossEntropyValidation"] = dist.crossEntropy(validationSamples)
        stats[i]["validationSamples"] = validationSamples
        stats[i]["samples"] = {"shuffled": {},
                               "not_shuffled": {}}
        stats[i]["samples"]["shuffled"]["rvs"] = dist.rvs(numSamples, shuffle=True)
        stats[i]["samples"]["shuffled"]["uniform_validation"] = dist.cdf(validationSamples,
                                                                         shuffle=True)
        kstests = [None] * numDims

        for idim in range(numDims):
            samples1d = stats[i]["samples"]["shuffled"]["uniform_validation"][:, idim]
            res_test = kstest(samples1d, Uniform(0, 1).cdf)
            kstests[idim] = res_test.statistic, res_test.pvalue
            if plot:
                plt.figure()
                plt.hist(samples1d, cumulative=True, normed=True)
                xs = np.linspace(0, 1, 10)
                plt.plot(xs, [Uniform(0, 1).cdf(xi) for xi in xs])
                plt.title("shuffled: %i, %s" % (idim, kstests[idim]))
        print("-" * 80)
        print("shuffled    ", kstests, np.min(kstests), np.max(kstests))
        if plot:
            plt.show()

        stats[i]["samples"]["shuffled"]["kstests"] = kstests
        stats[i]["samples"]["not_shuffled"]["rvs"] = dist.rvs(numSamples, shuffle=False)
        stats[i]["samples"]["not_shuffled"]["uniform_validation"] = dist.cdf(validationSamples,
                                                                             shuffle=False)
        kstests = [None] * numDims
        for idim in range(numDims):
            samples1d = stats[i]["samples"]["not_shuffled"]["uniform_validation"][:, idim]
            res_test = kstest(samples1d, Uniform(0, 1).cdf)
            kstests[idim] = res_test.statistic, res_test.pvalue
            if plot:
                plt.figure()
                plt.hist(samples1d, cumulative=True, normed=True)
                xs = np.linspace(0, 1, 1000)
                plt.plot(xs, [Uniform(0, 1).cdf(xi) for xi in xs])
                plt.title("not shuffled: %i, %s" % (idim, kstests[idim]))
        print("not shuffled", kstests, np.min(kstests), np.max(kstests))
        if plot:
            plt.show()

        stats[i]["samples"]["not_shuffled"]["kstests"] = kstests

        print("CV valid = %g" % crossEntropyValidation[i, 1])

        # write results to file
        if out:
            out_crossEntropy = os.path.join("data", method, "%s.%s.validation.cross_entropies.csv" % (method, functionName))
            np.savetxt(out_crossEntropy, crossEntropyValidation[:i, :])

            # save stats to pickle
            out_stats = os.path.join("data", method, "%s.%s.best.stats.pkl" % (method, functionName))
            fd = open(out_stats, "w")
            pkl.dump(stats, fd)
            fd.close()


def run_plotRoutine(tikz=False):
    stats = {}
    for config in density_configs:
        for root, _, files in os.walk(os.path.join("data", config)):
            for filename in files:
                if "pkl" in filename:
                    path = os.path.join(root, filename)
                    print("=" * 80)
                    print("plotting '%s'" % path)
                    print("=" * 80)
                    fd = open(path, "r")
                    currentStats = pkl.load(fd)
                    fd.close()
                    currentStats["config"]["functionName"] += "_i%i" % currentStats["config"]["iteration"]
                    currentStats["config"]["maxNumGridPoints"] = 0
                    if currentStats["config"]["candidateSearchAlgorithm"] == "intersections":
                        plotResults(currentStats, tikz)

                    stats[tuple(currentStats["config"].values())] = currentStats

        # plot joined stuff
        plotCostsPerIteration(stats)


# eval input parameters
density_configs = ["sgde_zero",
                   "sgde_boundaries",
                   "kde_gaussian",
                   "kde_epanechnikov",
                   "nataf"]

function_configs = ["mult_beta", "two_moons", "friedman"]


if __name__ == '__main__':
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--function', default="mult_beta", type=str, help='which data set should be used (mult_beta, mult_normal, mult_gamma, two_moons)')
    parser.add_argument('--method', default="sgde_zero", type=str, help='which method should be used to intepolate the new grid points (zero, log)')
    parser.add_argument('--candidates', default="join", type=str, help='which method should be used to intepolate the new grid points (zero, log)')
    parser.add_argument('--bandopt', default="rot", type=str, help='which method should be used to intepolate the new grid points (rot, ml)')
    parser.add_argument('--kfold', default=None, type=int, help='run kfold (20)')
    parser.add_argument('--numDims', default=2, type=int, help='number of dimensions for multivariate distributions')
    parser.add_argument('--numSamples', default=10000, type=int, help='number of samples that should be drawn from the best distributions')
    parser.add_argument('--plot', default=False, action='store_true', help='plot stuff')
    parser.add_argument('--out', default=False, action='store_true', help='write stuff to file')
    parser.add_argument('--tikz', default=True, action='store_true', help='save tikz images')
    args = parser.parse_args()

    if args.method not in density_configs:
        raise Exception("method '%s' is not known" % args.method)

    if args.kfold is not None:
        run_densityEstimation(args.function,
                              args.method,
                              kfold=args.kfold,
                              numDims=args.numDims,
                              numSamples=args.numSamples,
                              candidates=args.candidates,
                              bandwidthOptimizationType=args.bandopt,
                              out=args.out,
                              plot=args.plot,
                              tikz=args.tikz)
    else:
        run_plotRoutine(args.tikz)
