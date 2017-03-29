import os
import pickle as pkl
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.helper import sortPermutations
from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend, savefig


def get_key_sg(gridType, maxGridSize):
    return (gridType, maxGridSize)


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

                if currentStats["surrogate"] == "sg":
                    key = get_key_sg(currentStats["grid_type"],
                                     currentStats["max_grid_size"])
                    ans["sg"][key] = currentStats
                else:
                    key = get_key_mc(currentStats["sampling_strategy"],
                                     currentStats["num_model_evaluations"])
                    ans["mc"][key] = currentStats

                print "-" * 80
                print "loaded '%s'" % (key,)

    return ans


settings = {'beta': {'sg': [
                            ("polyBoundary", 10000),
                            ]}}

if __name__ == "__main__":
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="beta", type=str, help="define which probabilistic model should be used")
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    args = parser.parse_args()

    results = load_results(args.model)

    for error_type in ["l2test"]:
        # extract the ones needed for the table
        sg_settings = settings[args.model]["sg"]
        fig = plt.figure()
        if args.surrogate in ["sg", "both"]:
            for gridType, maxGridSize in sg_settings:
                key = get_key_sg(gridType, maxGridSize)
                n = len(results["sg"][key]["results"])
                num_evals = np.ndarray(n)
                errors = np.ndarray(n)
                for i, (boundaryLevel, values) in enumerate(results["sg"][key]["results"].items()):
                    num_evals[i] = boundaryLevel
                    errors[i] = values[error_type]
                ixs = np.argsort(num_evals)
                plt.plot(num_evals[ixs], errors[ixs], "o-",
                         label=r"SG (polyBoundary, $\ell=9$)")

        plt.xlim(9.5, 0.5)
        ticks = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.xticks(ticks, ticks)
        plt.yscale("log")
        plt.ylabel(r"$||u - u_{\mathcal{I}}||_{L_2(\Xi)}$")
        plt.xlabel(r"$\ell^{\text{b}}$")
        lgd = insert_legend(fig, loc="bottom")
        savefig(fig, "plots/sg_boundary_results")
