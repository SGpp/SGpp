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
    print "-" * 80
    print "load setting: %s" % setting

    filename = os.path.join("data", setting, "%s.best.stats.pkl" % setting)
    fd = open(filename, "r")
    stats = pkl.load(fd)
    fd.close()

    for iteration, values in stats.items():
        jsonObject = None
        if setting in ["sgde_zero", "sgde_boundaries"]:
            jsonObject = json.loads(values['posSGDE_json'])
        elif setting in ["kde_gaussian", "kde_epanechnikov"]:
            jsonObject = json.loads(values['KDEDist_json'])
        elif setting == "analytic_lognormalBeta":
            jsonObject = json.loads(values['Dist_json'])

        stats[iteration]["dist"] = Dist.fromJson(jsonObject)

    return stats


def plotLogLikelihood(densities):
    numDensities = len(densities)
    numIterations = len(densities.itervalues().next()) - 1

    data = np.ndarray((numIterations, numDensities))
    names = [None] * numDensities
    i = 0
    for i, (setting, stats) in enumerate(densities.items()):
        names[i] = setting.replace("_", "-")
        for j, values in enumerate(stats.values()):
            if j < 10:
                data[j, i] = values["crossEntropyValidation"]
            if j > 10:
                data[j - 1, i] = values["crossEntropyValidation"]

    pos = np.arange(0, numDensities)
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)
    plt.show()

def loadDensities():
    ans = {}
    for setting in density_configs:
        if "analytic" not in setting:
            ans[setting] = loadDensity(setting)
    return ans

if __name__ == "__main__":
    densities = loadDensities()
    plotLogLikelihood(densities)
