'''
Created on Jan 27, 2017

@author: franzefn
'''

import os
import pickle as pkl
import json
import numpy as np

import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.dists.Dist import Dist
from estimateDensity import density_configs


def loadDensity(setting):
    stats = None
    filename = os.path.join("data", setting, "%s.two_moons.best.stats.pkl" % setting)
    if os.path.exists(filename) and \
            "sgde" in setting:
        print "-" * 80
        print "load setting: %s" % setting
        fd = open(filename, "r")
        stats = pkl.load(fd)
        fd.close()

        for iteration, values in stats.items():
            jsonObject = None
            if setting in ["sgde_zero", "sgde_boundaries"]:
                jsonObject = json.loads(values['posSGDE_json'])
            elif setting in ["kde_gaussian", "kde_epanechnikov"]:
                jsonObject = json.loads(values['KDEDist_json'])
            elif setting in ["nataf"]:
                jsonObject = json.loads(values["NatafDist_json"])

            if jsonObject is not None:
                stats[iteration]["dist"] = Dist.fromJson(jsonObject)

    return stats


def plotLogLikelihood(densities):
    numDensities = len(densities)
    numIterations = len(densities.itervalues().next())

    data = np.ndarray((numIterations, numDensities))
    names = [None] * numDensities
    i = 0
    for i, (setting, stats) in enumerate(densities.items()):
        names[i] = setting.replace("_", "-")
        for j, values in enumerate(stats.values()):
            data[j, i] = values["crossEntropyValidation"]

    pos = np.arange(0, numDensities)
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)
    plt.show()

def plotpvalueofKolmogorovSmirnovTest(densities):
    numDensities = len(densities)
    numIterations = len(densities.itervalues().next())

    data = np.zeros((numIterations, 2 * numDensities))
    names = [None] * data.shape[1]
    i = 0
    for i, (setting, stats) in enumerate(densities.items()):
        names[2 * i] = setting.replace("_", "-") + " (shuffled)"
        names[2 * i + 1] = setting.replace("_", "-") + " (not shuffled)"
        for j, values in enumerate(stats.values()):
            numDims = values["config"]["numDims"]
            pvalues_shuffled = np.zeros(numDims)
            pvalues_not_shuffled = np.zeros(numDims)
            for idim in xrange(numDims):
                pvalues_shuffled[idim] = values["samples"]["shuffled"]["kstests"][idim][1]
                pvalues_not_shuffled[idim] = values["samples"]["not_shuffled"]["kstests"][idim][1]
            data[j, 2 * i] = pvalues_shuffled.mean()
            data[j, 2 * i + 1] = pvalues_not_shuffled.mean()

    pos = np.arange(0, len(names))
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)
    plt.ylabel(r"$p$-value of Kolmogorov-Smirnov test")
    plt.show()


def loadDensities():
    ans = {}
    for setting in density_configs:
        stats = loadDensity(setting)
        if stats is not None:
            ans[setting] = stats
    return ans

if __name__ == "__main__":
    densities = loadDensities()
#     plotLogLikelihood(densities)
    plotpvalueofKolmogorovSmirnovTest(densities)
