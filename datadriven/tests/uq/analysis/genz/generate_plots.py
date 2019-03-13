import os
import pickle as pkl
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.helper import sortPermutations
from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend, savefig, \
    load_font_properties, load_color, load_marker


def get_key_sg(gridType, maxGridSize, rosenblatt, boundaryLevel):
    return (gridType, maxGridSize, rosenblatt, boundaryLevel)


def get_key_mc(sampling_strategy, N):
    return (sampling_strategy, N)

def load_results(inputspace, path="results"):
    ans = {"sg": {},
           "mc": {}}
    for root, dirs, files in os.walk(os.path.join(path)):
        for filename in files:
            if "pkl" in filename:
                path = os.path.join(root, filename)
                fd = open(path, "r")
                currentStats = pkl.load(fd)
                fd.close()

                if "rosenblatt" not in currentStats:
                    currentStats["rosenblatt"] = False
                if "b3" in filename:
                    currentStats["boundaryLevel"] = 3
                if "b1" in filename:
                    currentStats["boundaryLevel"] = 1
                if "b10" in filename:
                    currentStats["boundaryLevel"] = 10

                if currentStats["surrogate"] == "sg":
                    key = get_key_sg(currentStats["grid_type"],
                                     currentStats["max_grid_size"],
                                     currentStats["rosenblatt"],
                                     currentStats["boundaryLevel"])
                    ans["sg"][key] = currentStats
                else:
                    key = get_key_mc(currentStats["sampling_strategy"],
                                     currentStats["num_model_evaluations"])
                    ans["mc"][key] = currentStats

                print("-" * 80)
                print("loaded '%s'" % (key,))

    return ans


def plotBoundaryResult(results):
    settings = {'beta': {'sg': [
                                ("polyBoundary", 10000, 10, "bound."),
                                ("polyClenshawCurtisBoundary", 10000, 10, "CC-bound.")
                                ]}}


    error_type = "l2test"
    # extract the ones needed for the table
    sg_settings = settings[args.model]["sg"]
    fig = plt.figure()
    for k, (gridType, maxGridSize, boundaryLevel, gridTypeLabel) in enumerate(sg_settings):
        key = get_key_sg(gridType, maxGridSize, False, boundaryLevel)
        n = len(results["sg"][key]["results"])
        num_evals = np.ndarray(n)
        errors = np.ndarray(n)
        for i, (boundaryLevel, values) in enumerate(results["sg"][key]["results"].items()):
            num_evals[i] = boundaryLevel
#             num_evals[i] = values["num_model_evaluations"]
            errors[i] = values[error_type]
        ixs = np.argsort(num_evals)
        plt.plot(num_evals[ixs], errors[ixs],
                 color=load_color(k),
                 marker=load_marker(k),
                 label=r"%s ($\ell=9$)" % gridTypeLabel)

    plt.xlim(9.5, 0.5)
    ticks = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    plt.xticks(ticks, ticks)
#     plt.xscale("log")
    plt.yscale("log")
    plt.ylabel(r"$||u - u_{\mathcal{I}}||_{L_2(\Xi)}$")
    plt.xlabel(r"$\ell^{\text{b}}$")
    plt.title(r"Regular SG (poly, $D=2$)",
              fontproperties=load_font_properties())
    lgd = insert_legend(fig, loc="bottom", ncol=1)
    savefig(fig, "plots/sg_boundary_results")


def plotConvergenceResults(results):
    settings = {'beta': {'sg': [
                                ("polyBoundary", 10000, False, 1, "bound."),
                                ("polyClenshawCurtisBoundary", 10000, False, 1, "CC-bound."),
#                                 ("modpoly", 10000, False, 1, "modified"),
#                                 ("modPolyClenshawCurtis", 10000, False, 1, "modified-CC"),
                                ("polyBoundary", 10000, True, 1, "bound."),
                                ("polyClenshawCurtisBoundary", 10000, True, 1, "CC-bound."),
                                ("poly", 10000, False, 1, "no bound.")
#                                 ("polyBoundary", 2000, False, 1, "poly-bound., exp"),
#                                 ("polyBoundary", 2001, False, 1, "poly-bound., l2"),
#                                 ("polyBoundary", 2002, False, 1, "poly-bound., simple")
                                ]}}


    error_type = "l2test"
    # extract the ones needed for the table
    sg_settings = settings[args.model]["sg"]
    fig = plt.figure()
    for k, (gridType, maxGridSize, rosenblatt, boundaryLevel, gridTypeLabel) in enumerate(sg_settings):
        key = get_key_sg(gridType, maxGridSize, rosenblatt, boundaryLevel)
        n = len(results["sg"][key]["results"])
        num_evals = np.ndarray(n)
        errors = np.ndarray(n)
        for i, (level, values) in enumerate(results["sg"][key]["results"].items()):
            num_evals[i] = values["num_model_evaluations"]
            errors[i] = values[error_type]
        print(num_evals)
        ixs = np.argsort(num_evals)
        if "bound" in gridTypeLabel and "no" not in gridTypeLabel:
            if rosenblatt:
                label = "%s ($\\ell^{\\text{b}}=%i$, Rosen.)" % (gridTypeLabel, boundaryLevel)
            else:
                label = r"%s ($\ell^{\text{b}}=%i$)" % (gridTypeLabel, boundaryLevel)
        else:
            if rosenblatt:
                label = r"%s (Rosenblatt)" % (gridTypeLabel,)
            else:
                label = r"%s" % (gridTypeLabel,)
        plt.loglog(num_evals[ixs], errors[ixs], "o-",
                   color=load_color(k),
                   marker=load_marker(k),
                   label=label)

    plt.ylabel(r"$||u - u_{\mathcal{I}}||_{L_2(\Xi)}$")
    plt.xlabel(r"\# number of grid points")
    plt.title(r"Regular SG (poly, $D=2$)",
              fontproperties=load_font_properties())
    lgd = insert_legend(fig, loc="bottom", ncol=1)
    savefig(fig, "plots/sg_convergence_results")

if __name__ == "__main__":
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="beta", type=str, help="define which probabilistic model should be used")
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    args = parser.parse_args()

    results = load_results(args.model)

    plotBoundaryResult(results)
    plotConvergenceResults(results)
