# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import os
import pickle as pkl
import json
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from estimateDensity import load_data_set
from scipy.stats import chi2_contingency

from pysgpp.extensions.datadriven.uq.dists.Dist import Dist
from estimateDensity import density_configs
from pysgpp.extensions.datadriven.uq.plot.colors import savefig, \
    load_font_properties, load_color, insert_legend
from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d


def loadDensity(setting, functionName):
    stats = None
    filename = os.path.join("data", setting, "%s.%s.best.stats.pkl" % (setting, functionName))
    if os.path.exists(filename) and \
            (("moons" in functionName and (("sgde" in setting) or
                                           "kde" in setting)) or
             ("beta" in functionName and (("sgde" in setting) or \
                                          "nataf" in setting or \
                                          "kde" in setting))):
        print("-" * 80)
        print("load setting: %s" % setting)
        fd = open(filename, "r")
        stats = pkl.load(fd)
        fd.close()

        for iteration, values in list(stats.items()):
            jsonObject = None
            if setting in ["sgde_zero", "sgde_boundaries"]:
                jsonObject = json.loads(values['posSGDE_json'])
            elif setting in ["kde_gaussian", "kde_epanechnikov"]:
                jsonObject = json.loads(values['KDEDist_json'])
            elif setting in ["nataf"]:
                jsonObject = json.loads(values["NatafDist_json"])
            if setting != "nataf" and jsonObject is not None:
                stats[iteration]["dist"] = Dist.fromJson(jsonObject)
        print("done")
    return stats


def plotLogLikelihood(densities, functionName, out=False):
    numDensities = len(densities)
    numIterations = 0
    for i, (setting, stats) in enumerate(densities.items()):
        numIterations = max(numIterations, len(stats))

    data = {"train": np.zeros((numIterations, numDensities)),
            "test": np.zeros((numIterations, numDensities)),
            "validation": np.zeros((numIterations, numDensities))}
    names = [None] * numDensities
    i = 0
    for i, setting in enumerate(["kde_gaussian",
                                 "kde_epanechnikov",
                                 "sgde_zero",
                                 "sgde_boundaries"]):
        stats = densities[setting]
        if "sgde" in setting:
            if "zero" in setting:
                names[i] = "SGDE \n set-to-zero"
            else:
                names[i] = "SGDE \n interp. bound."
            trainkey = "ZeroSGDE"
        elif "nataf" in setting:
            names[i] = "Nataf"
        elif "gaussian" in setting:
            names[i] = "KDE \n Gaussian"
            trainkey = "KDE"
        elif "epanechnikov" in setting:
            names[i] = "KDE \n Epan."
            trainkey = "KDE"
        for j, values in enumerate(stats.values()):
            data["train"][j, i] = values["crossEntropyTrain%s" % trainkey]
            data["test"][j, i] = values["crossEntropyTest%s" % trainkey]
            data["validation"][j, i] = values["crossEntropyValidation"]

    pos = np.arange(0, numDensities)

    fig = plt.figure(figsize=(13, 6.5))
    ax = fig.add_subplot(111)
#     plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
#                    showextrema=True, showmedians=True, bw_method=0.5)
    width = 0.28
    for i, category in enumerate(["train", "test", "validation"]):
        values = data[category]
        yval = np.ndarray(values.shape[1])
        for j in range(values.shape[1]):
            yval[j] = np.mean(values[:, j])
        rects = ax.bar(pos + i * width, yval, width, color=load_color(i),
                       label=category)
        for rect in rects:
            h = -rect.get_height()
            ax.text(rect.get_x() + (rect.get_width() / 2.), h - 0.2, '%.2f' % h,
                    ha='center', va='bottom')

#     plt.xticks(pos, names)
    ax.set_xticks(pos + width)
    ax.set_xticklabels(names)
    ax.set_ylabel("cross entropy")
    yticks = np.arange(-1.5, 0.5, 0.5)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set_ylim(-1.7, 0)
    lgd = insert_legend(fig, loc="right", ncol=1)

    if out:
        savefig(fig,
                os.path.join("plots", "log_likelihood_%s" % functionName),
                tikz=True,
                lgd=lgd)
        plt.close(fig)
    else:
        plt.show()


def plotpvalueofChi2IndependenceTest(densities,
                                     functionName,
                                     c=0.0,
                                     out=False):
    numDensities = len(densities)
    numIterations = 0
    for i, (setting, stats) in enumerate(densities.items()):
        numIterations = max(numIterations, len(stats))

    data = np.zeros((numIterations, 2 * numDensities))
    names = [None] * data.shape[1]
    i = 0
    for i, setting in enumerate(["kde_gaussian",
                                 "kde_epanechnikov",
                                 "sgde_zero",
                                 "sgde_boundaries"]):
        stats = densities[setting]
        if "sgde" in setting:
            if "zero" in setting:
                names[2 * i] = "SGDE \n set-to-zero \n shuffled"
                names[2 * i + 1] = "SGDE \n set-to-zero \n not shuffled"
            else:
                names[2 * i] = "SGDE \n interp. bound. \n shuffled"
                names[2 * i + 1] = "SGDE \n interp. bound. \n not shuffled"
        elif "nataf" in setting:
            names[2 * i] = "Nataf \n shuffled"
            names[2 * i + 1] = "Nataf \n not shuffled"
        elif "gaussian" in setting:
            names[2 * i] = "KDE \n Gaussian \n shuffled"
            names[2 * i + 1] = "KDE \n Gaussian \n not shuffled"
        elif "epanechnikov" in setting:
            names[2 * i] = "KDE \n Epan. \n shuffled"
            names[2 * i + 1] = "KDE \n Epan. \n not shuffled"
        for j, values in enumerate(stats.values()):
            numDims = values["config"]["numDims"]

            # apply the chi 2 test
            bins = np.linspace(0, 1, 10)
            samples = values["samples"]["shuffled"]["uniform_validation"]
            inner_samples = np.array([])
            for sample in samples:
                if c < sample[0] < 1 - c and c < sample[1] < 1 - c:
                    inner_samples = np.append(inner_samples, sample)
            inner_samples = inner_samples.reshape((inner_samples.size // 2), 2)
            h0 = np.histogram2d(inner_samples[:, 0],
                                inner_samples[:, 1],
                                bins=bins)[0][2:-2, 2:-2]
            pvalue_shuffled = chi2_contingency(h0)[1]

            if False and j == 0:
                plt.figure()
                plt.scatter(inner_samples[:, 0], inner_samples[:, 1])

                plt.figure()
                plt.hist2d(inner_samples[:, 0], inner_samples[:, 1], bins=20)
                plt.colorbar()
                plt.title("%s shuffled, %g" % (setting.replace("_", " "),
                                               pvalue_shuffled))

            samples = values["samples"]["not_shuffled"]["uniform_validation"]
            inner_samples = np.array([])
            for sample in samples:
                if c < sample[0] < 1 - c and c < sample[1] < 1 - c:
                    inner_samples = np.append(inner_samples, sample)
            inner_samples = inner_samples.reshape((inner_samples.size // 2), 2)
            h0 = np.histogram2d(inner_samples[:, 0],
                                inner_samples[:, 1],
                                bins=bins)[0][2:-2, 2:-2]
            pvalue_not_shuffled = chi2_contingency(h0)[1]

            if False and j == 0:
                plt.figure()
                plt.scatter(inner_samples[:, 0],
                            inner_samples[:, 1])

                plt.figure()
                plt.hist2d(inner_samples[:, 0],
                           inner_samples[:, 1], bins=20)
                plt.colorbar()
                plt.title("%s not shuffled, %g" % (setting.replace("_", " "),
                                                   pvalue_not_shuffled))

                plt.show()

            data[j, 2 * i] = pvalue_shuffled
            data[j, 2 * i + 1] = pvalue_not_shuffled

    pos = np.arange(0, len(names))
    xlim = (np.min(pos) - 0.5, np.max(pos) + 0.5)
    fig = plt.figure(figsize=(17, 5))
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)
    plt.ylabel("$p$-value")
    plt.hlines(0.05, xlim[0], xlim[1], linestyle="--")
    plt.xlim(xlim)

    if "moons" in functionName:
        plt.title("$\chi^2$ test",
                  fontproperties=load_font_properties())
    else:
        plt.title("$\chi^2$ test",
                  fontproperties=load_font_properties())

    if out:
        savefig(fig,
                os.path.join("plots", "chi_squared_%s_c%i" % (functionName, 
                                                              np.round(c * 100))),
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
    for i, setting in enumerate(["kde_gaussian",
                                 "kde_epanechnikov",
                                 "sgde_zero",
                                 "sgde_boundaries"]):
        stats = densities[setting]
        if "sgde" in setting:
            if "zero" in setting:
                names[2 * i] = "SGDE \n set-to-zero \n shuffled"
                names[2 * i + 1] = "SGDE \n set-to-zero \n not shuffled"
            else:
                names[2 * i] = "SGDE \n interp. bound. \n shuffled"
                names[2 * i + 1] = "SGDE \n interp. bound. \n not shuffled"
        elif "nataf" in setting:
            names[2 * i] = "Nataf \n shuffled"
            names[2 * i + 1] = "Nataf \n not shuffled"
        elif "gaussian" in setting:
            names[2 * i] = "KDE \n Gaussian \n shuffled"
            names[2 * i + 1] = "KDE \n Gaussian \n not shuffled"
        elif "epanechnikov" in setting:
            names[2 * i] = "KDE \n Epan. \n shuffled"
            names[2 * i + 1] = "KDE \n Epan. \n not shuffled"
        for j, values in enumerate(stats.values()):
            numDims = values["config"]["numDims"]
            pvalues_shuffled = np.zeros(numDims)
            pvalues_not_shuffled = np.zeros(numDims)
            for idim in range(numDims):
                pvalues_shuffled[idim] = values["samples"]["shuffled"]["kstests"][idim][1]
                pvalues_not_shuffled[idim] = values["samples"]["not_shuffled"]["kstests"][idim][1]
            data[j, 2 * i] = pvalues_shuffled.mean()
            data[j, 2 * i + 1] = pvalues_not_shuffled.mean()

    pos = np.arange(0, len(names))
    xlim = (np.min(pos) - 0.5, np.max(pos) + 0.5)
    fig = plt.figure(figsize=(17, 5))
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)
    plt.ylabel("$p$-value")
    plt.hlines(0.05, xlim[0], xlim[1], linestyle="--")
    plt.xlim(xlim)


    if "moons" in functionName:
        plt.title("Kolmogorov-Smirnov test",
                  fontproperties=load_font_properties())
    else:
        plt.title("Kolmogorov-Smirnov test",
                  fontproperties=load_font_properties())

    if out:
        savefig(fig,
                os.path.join("plots", "kolmogorov_smirnov_%s" % functionName),
                tikz=True)
        plt.close(fig)
    else:
        plt.show()

def plotCovarianceConvergence(densities, functionName, out=False):
    _, _, natafType = load_data_set(functionName, numSamples=0, numDims=2)

    covMatrix = natafType["cov"]

    numDensities = len(densities)
    numIterations = 0
    for i, (setting, stats) in enumerate(densities.items()):
        numIterations = max(numIterations, len(stats))

    data = np.zeros((numIterations, numDensities))
    names = [None] * numDensities
    i = 0
    for i, setting in enumerate(["kde_gaussian",
                                 "kde_epanechnikov",
                                 "sgde_zero",
                                 "sgde_boundaries"]):
        stats = densities[setting]
        if "sgde" in setting:
            if "zero" in setting:
                names[i] = "SGDE \n set-to-zero"
            else:
                names[i] = "SGDE \n interp. bound."
        elif "nataf" in setting:
            names[i] = "Nataf"
        elif "gaussian" in setting:
            names[i] = "KDE \n Gaussian"
        elif "epanechnikov" in setting:
            names[i] = "KDE \n Epan."
        for j, values in enumerate(stats.values()):
            data[j, i] = np.linalg.norm(values["dist"].cov() - covMatrix)

    pos = np.arange(0, numDensities)

    fig = plt.figure(figsize=(13, 4.5))
    ax = fig.add_subplot(111)
    plt.violinplot(data, pos, points=60, widths=0.7, showmeans=True,
                   showextrema=True, showmedians=True, bw_method=0.5)
    plt.xticks(pos, names)
    plt.ylabel(r"$||\hat{C} - C||$")

    if out:
        savefig(fig,
                os.path.join("plots", "convergence_covariance_%s" % functionName),
                tikz=True)
        plt.close(fig)
    else:
        plt.show()


def plotDataset(functionName, numSamples=10000, numDims=2, out=False):
    dataset, bounds, _ = load_data_set(functionName, numSamples, numDims=2)
    fig = plt.figure()
    plt.plot(dataset[:, 0], dataset[:, 1],
             "o ",
             color=load_color(0))
    plt.xlabel(r"$\xi_1$")
    plt.ylabel(r"$\xi_2$")
    plt.xlim(bounds[0])
    plt.ylim(bounds[1])
    xticks = np.arange(0, 1.2, 0.2)
    plt.xticks(xticks, [str(xi) for xi in xticks])
    plt.yticks(xticks, [str(xi) for xi in xticks])
    plt.title("Two-moons dataset",
              fontproperties=load_font_properties())

    if out:
        filename = os.path.join("plots", "%s_dataset" % functionName)
        print(filename)
        fig.set_size_inches(5.7, 5, forward=True)
        savefig(fig, filename, tikz=True)
        plt.close(fig)
    else:
        plt.show()

def plotDensities(densities, functionName, out=False):
    for setting, stats in list(densities.items()):
        if setting != "nataf":
            U = stats[0]["dist"]
            if "kde" in setting:
                label = r'$f_{\mathcal{S}_M}^{\kappa}(\xi_1, \xi_2)$'
                if "gaussian" in setting:
                    title = "KDE (Gaussian)"
                else:
                    title = "KDE (Epanechnikov)"
            else:
                label = r'$f_{\mathcal{I}}(\xi_1, \xi_2)$'
                if "zero" in setting:
                    title = "SGDE (set-to-zero)"
                else:
                    title = "SGDE (interp. bound.)"
            fig = plt.figure()
            plotDensity2d(U, color_bar_label=label)
            plt.xlabel(r"$\xi_1$")
            plt.ylabel(r"$\xi_2$")
            xticks = np.arange(0, 1.2, 0.2)
            plt.xticks(xticks, [str(xi) for xi in xticks])
            plt.yticks(xticks, [str(xi) for xi in xticks])
            plt.title(title,
                      fontproperties=load_font_properties())
            if out:
                filename = os.path.join("plots", "%s_%s" % (functionName, setting))
                print(filename)
                fig.set_size_inches(5.7, 5, forward=True)
                savefig(fig, filename, tikz=True)
                plt.close(fig)
            
    if not out:
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
    plotpvalueofChi2IndependenceTest(densities, args.function, c=0.0, out=args.out)
#     plotpvalueofChi2IndependenceTest(densities, args.function, c=0.25, out=args.out)
#     plotLogLikelihood(densities, args.function, args.out)
    plotpvalueofKolmogorovSmirnovTest(densities, args.function, args.out)
#     plotDensities(densities, args.function, out=args.out)
#     plotDataset(args.function, out=args.out)
#     plotCovarianceConvergence(densities, "mult_beta", args.out)
