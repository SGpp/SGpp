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
from pysgpp.extensions.datadriven.uq.plot.colors import savefig, \
    load_font_properties
from argparse import ArgumentParser


def loadDensity(setting, functionName):
    stats = None
    filename = os.path.join("data", setting, "%s.%s.best.stats.pkl" % (setting, functionName))
    if os.path.exists(filename) and \
            (("moons" in functionName and (("sgde" in setting and "zero" in setting) or
                                           "kde" in setting and "gaussian" in setting)) or
             ("beta" in functionName and (("sgde" in setting and "zero" in setting) or \
                                          "nataf" in setting or \
                                          "kde" in setting and "gaussian" in setting))):
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


def plotLogLikelihood(densities, functionName, out=False):
    numDensities = len(densities)
    numIterations = 0
    for i, (setting, stats) in enumerate(densities.items()):
        numIterations = max(numIterations, len(stats))

    data = np.zeros((numIterations, numDensities))
    names = [None] * numDensities
    i = 0
    for i, (setting, stats) in enumerate(densities.items()):
        if "sgde" in setting:
            if "zero" in setting:
                names[i] = "SGDE \n extended"
        elif "nataf" in setting:
            names[i] = "Nataf"
        elif "gaussian" in setting:
            names[i] = "KDE \n Gaussian"
        for j, values in enumerate(stats.values()):
            data[j, i] = values["crossEntropyValidation"]

    pos = np.arange(0, numDensities)

    fig = plt.figure()
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)

    if out:
        savefig(fig,
                os.path.join("plots", "log_likelihood_%s" % functionName),
                tikz=True)
        plt.close(fig)
    else:
        plt.show()


def plotpvalueofKolmogorovSmirnovTest(densities, functionName, out=False):
    numDensities = len(densities)
    numIterations = 0
    for i, (setting, stats) in enumerate(densities.items()):
        numIterations = max(numIterations, len(stats))

    data = np.zeros((numIterations, 2 * numDensities))
    names = [None] * data.shape[1]
    i = 0
    for i, (setting, stats) in enumerate(densities.items()):
        if "sgde" in setting:
            if "zero" in setting:
                names[2 * i] = "SGDE \n extended \n shuffled"
                names[2 * i + 1] = "SGDE \n extended \n not shuffled"
            else:
                names[2 * i] = "SGDE \n interp. bound. \n shuffled"
                names[2 * i + 1] = "SGDE \n interp. bound. \n not shuffled"
        elif "nataf" in setting:
            names[2 * i] = "Nataf \n shuffled"
            names[2 * i + 1] = "Nataf \n not shuffled"
        elif "gaussian" in setting:
            names[2 * i] = "KDE \n Gaussian \n shuffled"
            names[2 * i + 1] = "KDE \n Gaussian \n not shuffled"
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
    fig = plt.figure()
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)
    plt.ylabel("$p$-value")
    plt.title("Two moons: \n Kolmogorov-Smirnov test",
              fontproperties=load_font_properties())

    if out:
        savefig(fig,
                os.path.join("plots", "kolmogorov_smirnov_%s" % functionName),
                tikz=True)
        plt.close(fig)
    else:
        plt.show()


def loadDensities(functionName):
    ans = {}
    for setting in density_configs:
        stats = loadDensity(setting, functionName)
        if stats is not None:
            ans[setting] = stats
    return ans

if __name__ == "__main__":
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--function', type=str, default="two_moons", help='function to be interpolated (normal, beta)')
    parser.add_argument('--out', default=False, action='store_true', help='write stuff to file')
    args = parser.parse_args()

    densities = loadDensities(args.function)
    plotLogLikelihood(densities, args.function, args.out)
    plotpvalueofKolmogorovSmirnovTest(densities, args.function, args.out)
    plt.show()
