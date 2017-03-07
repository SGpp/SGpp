'''
Created on Apr 18, 2016

@author: franzefn
'''

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
    createOperationMakePositive, DataVector
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d

from positivity.plot import plotResults, plotCostsPerIteration
from argparse import ArgumentParser
from positivity.helper_functions import strToCandidateSearchAlgorithm, \
    strToInterpolationAlgorithm


def splitset(samples, splitPercentage=0.8):
    numSamples = samples.shape[0]
    numTrainSamples = int(np.ceil(numSamples * splitPercentage))
    trainSamplesixs = np.random.choice(numSamples, numTrainSamples, replace=False)

    trainSamplesMask = np.zeros(numSamples, dtype="bool")
    trainSamplesMask[trainSamplesixs] = True
    testSamplesMask = np.invert(trainSamplesMask)

    assert numTrainSamples == samples[trainSamplesMask].shape[0]
    
    return samples[trainSamplesMask], samples[testSamplesMask]


def estimateSGDEDensity(trainSamples, testSamples=None, iteration=0, plot=False, out=True, label="sgde_zero",
                        candidates="intersections", interpolation="setToZero"):
    print "train: %i x %i (mean=%g, var=%g)" % (trainSamples.shape[0], trainSamples.shape[1], np.mean(trainSamples), np.var(trainSamples))
    if testSamples is not None:
        print "test : %i x %i (mean=%g, var=%g)" % (testSamples.shape[0], testSamples.shape[1], np.mean(testSamples), np.var(testSamples))

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
    stats = {'config': {'functionName': 'twoMoons',
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

    for level in xrange(2, 5):
        print "-" * 60
        print "l=%i" % level
        for refinementSteps in xrange(0, 6):
            config["grid_level"] = level
            config["refinement_numSteps"] = refinementSteps
            sgdeDist = SGDEdist.byLearnerSGDEConfig(trainSamples, config=config)
            # -----------------------------------------------------------
            grid, alpha = sgdeDist.grid, sgdeDist.alpha
            cvSgde = sgdeDist.crossEntropy(testSamples)

            maxLevel = grid.getStorage().getMaxLevel()
            numDims = grid.getStorage().getDimension()

            print "  " + "-" * 30
            print "  #ref = %i: gs=%i -> CV test = %g" % (refinementSteps,
                                                          sgdeDist.grid.getSize(),
                                                          cvSgde)
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
            positiveSgdeDist = SGDEdist(positiveGrid, positiveAlpha, trainSamples)
            # -----------------------------------------------------------
            cvPositiveSgde = positiveSgdeDist.crossEntropy(testSamples)

            if plot:
                fig = plt.figure()
                plotSG2d(grid, alpha, show_negative=True, show_grid_points=True)
                plt.title("pos: N=%i: vol=%g, log=%g" % (positiveGrid.getSize(),
                                                         doQuadrature(positiveGrid, positiveAlpha),
                                                         cvPositiveSgde))
                fig.show()

            # -----------------------------------------------------------
            print "  positive: gs=%i -> CV test = %g" % (positiveGrid.getSize(), cvPositiveSgde)
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
                M = np.sum([1 for i in xrange(len(alpha)) if alpha[i] < 0])
                for d in xrange(2, numDims + 1):
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

                if plot:
                    if numDims == 2:
                        fig = plt.figure()
                        plotSG2d(positiveGrid, positiveAlpha, show_negative=True, show_grid_points=False,
                                 colorbarLabel=r"$\hat{f}_{\mathcal{I}_\text{SG}, \mathcal{I}_\text{ext}}$")
                        plt.title(r"positive: $N=%i/%i$; \# comparisons$=%i$" % (positiveGrid.getSize(),
                                                                                 (2 ** maxLevel - 1) ** numDims,
                                                                                 numComparisons))
                        plt.xlabel(r"porosity $\phi$")
                        plt.ylabel(r"permeability $\log(K_a)$")
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
                        ax.set_zlabel(r"$\hat{f}_{\mathcal{I}_\text{SG}, \mathcal{I}_\text{ext}}$", fontsize=20)
                        ax.set_xlabel(r"$\xi_1$", fontsize=20)
                        ax.set_ylabel(r"$\xi_2$", fontsize=20)

                        plt.tight_layout()
                        plt.savefig(os.path.join(pathResults, "%s_pos_i%i_l%i_r%i_3d.jpg" % (label, iteration, level, refinementSteps)))
                        plt.savefig(os.path.join(pathResults, "%s_pos_i%i_l%i_r%i_3d.pdf" % (label, iteration, level, refinementSteps)))
                        if out:
                            plt.close(fig)

            if plot:
                plt.show()


    if out:
        # save stats
        filename = os.path.join("data", label, "stats_d%i_a%i_r%i_i%i_%s_%s.pkl" % (numDims, 1, refinementSteps, iteration, candidates, interpolation))
        fd = open(filename, "w")
        pkl.dump(stats, fd)
        fd.close()
        print "stats saved to -> '%s'" % filename

        # dictionary that stores the information on the estimated densities
        myjson = {"Grid": {"dimNames": ["phi", "log(K_A)"],
                           "matrixEntries": ["phi", "log(K_A)"]},
                  "Set": {"path": "",
                          "grids": [],
                          "alphas": [],
                          "paramValues": [],
                          "paramName": "grid_size"}}
    
        for key, result in results.items():
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
        for out in crossEntropies.values():
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


def estimateKDEDensity(trainSamples, testSamples=None, iteration=0, plot=False, out=True, label="kde_gaussian"):
    print "train: %i x %i (mean=%g, var=%g)" % (trainSamples.shape[0], trainSamples.shape[1], np.mean(trainSamples), np.var(trainSamples))
    if testSamples is not None:
        print "test : %i x %i (mean=%g, var=%g)" % (testSamples.shape[0], testSamples.shape[1], np.mean(testSamples), np.var(testSamples))

    if "gaussian" in label:
        kernelType = KernelType_GAUSSIAN
    elif "epanechnikov" in label:
        kernelType = KernelType_EPANECHNIKOV
    else:
        raise AttributeError("label is unknown")

    kdeDist = KDEDist(trainSamples,
                      kernelType=kernelType,
                      bandwidthOptimizationType=BandwidthOptimizationType_MAXIMUMLIKELIHOOD)
    # -----------------------------------------------------------
    cvKDE = kdeDist.crossEntropy(testSamples)

    if plot:
        fig = plt.figure()
        plotDensity2d(kdeDist)
        plt.title("log=%g" % cvKDE)
        fig.show()
        plt.show()

    print "CV test = %g" % cvKDE

    # -----------------------------------------------------------
    if out:        
        pathResults = os.path.join("data", label)
    
        # serialize cross entropies
        out_crossEntropies = os.path.join(pathResults, "kde_cross_entropies.i%i.csv" % iteration)
        fd = open(out_crossEntropies, 'wb')
        file_writer = csv.writer(fd)
        file_writer.writerow(["crossEntropy"])
        file_writer.writerow([cvKDE])
        fd.close()
    
        # serialize samples
        np.savetxt(os.path.join(pathResults, "kde_train_samples.i%i.csv" % iteration), trainSamples)
        np.savetxt(os.path.join(pathResults, "kde_test_samples.i%i.csv" % iteration), testSamples)

        if plot:
            # plot density
            fig = plt.figure()
            plotDensity2d(kdeDist)
            plt.title("%s -> CV = %g" % (kdeDist.getBandwidths(), cvKDE))
            plt.savefig(os.path.join(pathResults, "kde_pdf.i%i.jpg" % iteration))
            plt.close(fig)

        # serialize best configuration to json
        out_bestDist = os.path.join(pathResults, "kde_best_config.i%i.json" % iteration)
        text = kdeDist.toJson()
        fd = open(out_bestDist, "w")
        fd.write(text)
        fd.close()

    # stats
    stats = {'config': {'functionName': 'twoMoons',
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


def run_densityEstimation(method,
                          kfold=20,
                          numSamples=1000,
                          candidates="join",
                          out=True,
                          tikz=False):
    if method == "sgde_zero":
        interpolation = "zero"
    else:  # interpolation == "boundaries":
        interpolation = "boundaries"

    samples = np.loadtxt("data/twomoons.csv")

    # do kfold cross validation
    crossEntropyValidation = np.zeros((kfold, 2))
    np.random.seed(1234567)
    learnSamples, validationSamples = splitset(samples, splitPercentage=0.7)

    stats = {}
    for i in xrange(kfold):
        print "=" * 100
        print "run (%s)= %i/%i" % (method, i + 1, kfold)
        print "=" * 100
        print "valid: %i x %i (mean=%g, var=%g)" % (validationSamples.shape[0], validationSamples.shape[1], np.mean(validationSamples), np.var(validationSamples))

        np.random.seed(i * 123456 + i % 2)
        trainSamples, testSamples = splitset(learnSamples, splitPercentage=1. - 1. / kfold)

        if "sgde" in method:
            dist, stats[i] = estimateSGDEDensity(trainSamples, testSamples, iteration=i, plot=plot, label=method, out=out,
                                                 candidates=candidates, interpolation=interpolation)
        elif "kde" in method:
            dist, stats[i] = estimateKDEDensity(trainSamples, testSamples, iteration=i, plot=plot, label=method, out=out)
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
        stats[i]["samples"]["shuffled"]["uniform_validation"] = dist.cdf(validationSamples, shuffle=True)
        stats[i]["samples"]["not_shuffled"]["rvs"] = dist.rvs(numSamples, shuffle=False)
        stats[i]["samples"]["not_shuffled"]["uniform_validation"] = dist.cdf(validationSamples, shuffle=False)

        print "CV valid = %g" % crossEntropyValidation[i, 1]

        # write results to file
        if out:
            out_crossEntropy = os.path.join("data", method, "%s.validation.cross_entropies.csv" % method)
            np.savetxt(out_crossEntropy, crossEntropyValidation[:i, :])

            # save stats to pickle
            out_stats = os.path.join("data", method, "%s.best.stats.pkl" % method)
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
                    print "=" * 80
                    print "plotting '%s'" % path
                    print "=" * 80
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
                   "kde_epanechnikov"]

if __name__ == '__main__':
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--method', default="sgde_zero", type=str, help='which method should be used to intepolate the new grid points (zero, log)')
    parser.add_argument('--candidates', default="join", type=str, help='which method should be used to intepolate the new grid points (zero, log)')
    parser.add_argument('--kfold', default=None, type=int, help='run kfold (20)')
    parser.add_argument('--numSamples', default=10000, type=int, help='number of samples that should be drawn from the best distributions')
    parser.add_argument('--plot', default=False, action='store_true', help='plot stuff')
    parser.add_argument('--out', default=False, action='store_true', help='write stuff to file')
    parser.add_argument('--tikz', default=True, action='store_true', help='save tikz images')
    args = parser.parse_args()

    if args.method not in density_configs:
        raise Exception("method '%s' is not known" % args.method)

    if args.kfold is not None:
        run_densityEstimation(args.method,
                              kfold=args.kfold,
                              numSamples=args.numSamples,
                              candidates=args.candidates,
                              out=args.out,
                              plot=args.plot,
                              tikz=args.tikz)
    else:
        run_plotRoutine(args.tikz)
