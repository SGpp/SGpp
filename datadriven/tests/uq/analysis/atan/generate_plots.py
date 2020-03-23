# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import os
import pickle as pkl
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.helper import sortPermutations
from argparse import ArgumentParser


def get_key_pce(expansion, sampling_strategy, N):
    return (expansion, sampling_strategy, N)

def get_key_sg(gridType, level, maxGridSize, refinement, isFull):
    return (gridType, level, maxGridSize, refinement, isFull)


def get_key_mc(sampling_strategy, N):
    return (sampling_strategy, N)

def load_results(inputspace, path="results"):
    ans = {"pce": {},
           "sg": {},
           "mc": {}}
    for root, dirs, files in os.walk(os.path.join(path, inputspace)):
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

                print("-" * 80)
                print("loaded '%s'" % (key,))

    return ans


settings = {'uniform': {'sg': [
#                                ("polyBoundary", 0, 3000, False, False),
#                                ("linearBoundary", 0, 3000, False, False),
#                                ("modLinearClenshawCurtis", 0, 3000, False, False),
#                                ("modlinear", 0, 3000, False, False),
#                                ("modPolyClenshawCurtis", 0, 3000, False, False),
#                                ("modpoly", 0, 3000, False, False),
                                ("linearClenshawCurtisBoundary", 0, 3000, False, False),
                                ("polyClenshawCurtisBoundary", 0, 3000, False, False)
#                                # ---------------------------------------------
#                                ("linearBoundary", 2, 3000, "var", False),
#                                ("polyBoundary", 2, 3000, "var", False),
#                                ("polyBoundary", 2, 3000, "squared", False),
#                                ("polyBoundary", 2, 3000, "weighted", False),
#                                ("polyBoundary", 2, 3000, "simple", False),
#                                ("polyBoundary", 2, 3000, "exp", False),
                               ],
                        'pce': [("full_tensor", 'gauss', 4000),
                                ('total_degree', 'gauss_leja', 4000)]},
            'beta': {'sg': [("polyBoundary", 0, 3000, False, False),
                            ("linearBoundary", 0, 3000, False, False),
                            # ---------------------------------------------
                            ("linearBoundary", 2, 3000, "var", False),
                            ("polyBoundary", 2, 3000, "var", False),
                            ("polyBoundary", 2, 3000, "squared", False),
                            ("polyBoundary", 2, 3000, "weighted", False),
                            ("polyBoundary", 2, 3000, "simple", False),
                            ("polyBoundary", 2, 3000, "exp", False)
                            ],
                     'pce': [("full_tensor", 'gauss', 4000),
                             ('total_degree', 'gauss_leja', 4000)]}}

if __name__ == "__main__":
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="uniform", type=str, help="define which probabilistic model should be used")
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    args = parser.parse_args()

    results = load_results(args.model)

    for error_type in ["l2test", "mean_error", "var_error"]:
        # extract the ones needed for the table
        pce_settings = settings[args.model]["pce"]
        sg_settings = settings[args.model]["sg"]
        plt.figure()
        if args.surrogate in ["pce", "both"]:
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

        if args.surrogate in ["sg", "both"]:
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

        plt.title(error_type.replace("_", " "))
        plt.legend(loc="lower left")
        plt.show()
