import os
import pickle as pkl
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.helper import sortPermutations


def get_key_pce(expansion, sampling_strategy, N):
    return (expansion, sampling_strategy, N)

def get_key_sg(gridType, level, maxGridSize, refinement, isFull):
    return (gridType, level, maxGridSize, refinement, isFull)


def get_key_mc(sampling_strategy, N):
    return (sampling_strategy, N)

def load_results(path="results"):
    ans = {"pce": {},
           "sg": {},
           "mc": {}}
    for root, dirs, files in os.walk(path):
        for filename in files:
            if "pkl" in filename:
                path = os.path.join(root, filename)
                fd = open(path, "r")
                currentStats = pkl.load(fd)
                fd.close()

                if currentStats["surrogate"] == "pce":
                    key = get_key_pce(currentStats["expansion"],
                                      currentStats["sampling_strategy"],
                                      currentStats["max_num_samples"])
                    ans["pce"][key] = currentStats
                elif currentStats["surrogate"] == "sg":
                    level = 0
                    if currentStats["refinement"]:
                        level = currentStats["level"]
                    key = get_key_sg(currentStats["grid_type"],
                                     level,
                                     currentStats["max_grid_size"],
                                     currentStats["refinement"],
                                     currentStats["is_full"])
                    ans["sg"][key] = currentStats
                else:
                    key = get_key_mc(currentStats["sampling_strategy"],
                                     currentStats["num_model_evaluations"])
                    ans["mc"][key] = currentStats

                print "-" * 80
                print "loaded '%s'" % (key,)

    return ans


if __name__ == "__main__":
    results = load_results()

    for error_type in ["l2test", "mean_error", "var_error"]:
        # extract the ones needed for the table
        pce_settings = [("full_tensor", 'gauss', 4000),
                        ('total_degree', 'gauss_leja', 4000)]
        sg_settings = [("polyBoundary", 0, 3000, False, True),
                       ("polyBoundary", 0, 3000, False, False),
                       ("polyBoundary", 2, 3000, "var", False),
                       ("polyBoundary", 2, 3000, "squared", False),
                       ("polyBoundary", 2, 3000, "weighted", False),
                       ("polyBoundary", 2, 3000, "exp", False)]

        plt.figure()
        for expansion, sampling_strategy, N in pce_settings:
            key = get_key_pce(expansion, sampling_strategy, N)
            n = len(results["pce"][key]["results"])
            num_evals = np.ndarray(n)
            errors = np.ndarray(n)
            for i, (num_samples, values) in enumerate(results["pce"][key]["results"].items()):
                num_evals[i] = num_samples
                errors[i] = values[error_type]
            ixs = np.argsort(num_evals)
            plt.loglog(num_evals[ixs], errors[ixs], "o-",
                       label=("pce (%s, %s)" % (expansion, sampling_strategy)).replace("_", " "))

        for gridType, level, maxGridSize, refinement, isFull in sg_settings:
            key = get_key_sg(gridType, level, maxGridSize, refinement, isFull)
            n = len(results["sg"][key]["results"])
            num_evals = np.ndarray(n)
            errors = np.ndarray(n)
            for i, values in enumerate(results["sg"][key]["results"].values()):
                num_evals[i] = values["num_model_evaluations"]
                errors[i] = values[error_type]
            ixs = np.argsort(num_evals)
            plt.loglog(num_evals[ixs], errors[ixs], "o-",
                       label=("sg (%s, %s)" % (gridType, refinement)).replace("_", " "))

        plt.title(error_type)
        plt.legend(loc="lower left")
        plt.show()
