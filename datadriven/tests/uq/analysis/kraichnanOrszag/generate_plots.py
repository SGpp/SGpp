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

def load_results(inputspace, setting, qoi, path="results"):
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

                if currentStats["qoi"] == qoi and currentStats["setting"] == setting:
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


settings = {'uniform': {'mc': [('latin_hypercube', 1000)],
                        'sg': [(8, 0, 1025, None, False)],
                        'pce': [("full_tensor", 'gauss', 4000),
                                ('total_degree', 'gauss_leja', 4000)]},
            'beta': {'sg': [
                            ("polyBoundary", 0, 3000, False, False)
                            ],
                     'pce': [("full_tensor", 'gauss', 4000),
                             ('total_degree', 'gauss_leja', 4000)]}}

if __name__ == "__main__":
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="uniform", type=str, help="define which probabilistic model should be used")
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--setting', default=1, type=int, help='parameter settign for test problem')
    parser.add_argument('--qoi', default="y2", type=str, help="define the quantity of interest")
    args = parser.parse_args()

    results = load_results(args.model, args.setting, args.qoi)

    for error_type in ["var"]:
        # extract the ones needed for the table
        mc_settings = settings[args.model]["mc"]
        pce_settings = settings[args.model]["pce"]
        sg_settings = settings[args.model]["sg"]
        plt.figure()

        # plot mc results to compare
        for sampling_strategy, numSamples in mc_settings:
            key = get_key_mc(sampling_strategy, numSamples)
            res = {'t': np.array([]),
                   'mean': np.array([]),
                   'var': np.array([])}
            for t, values in results["mc"][key]["results"].items():
                res['t'] = np.append(res["t"], t)
                res['mean'] = np.append(res["mean"], values["mean_estimated"])
                res['var'] = np.append(res["var"], values["var_estimated"])

            ixs = np.argsort(res["t"])
            plt.plot(res['t'][ixs], res[error_type][ixs], "-",
                     linewidth=2,
                     label="mc (N=%i)" % (numSamples,))

#         if args.surrogate in ["pce", "both"]:
#             for expansion, sampling_strategy, N in pce_settings:
#                 key = get_key_pce(expansion, sampling_strategy, N)
#                 n = len(results["pce"][key]["results"])
#                 num_evals = np.ndarray(n)
#                 errors = np.ndarray(n)
#                 for i, (num_samples, values) in enumerate(results["pce"][key]["results"].items()):
#                     num_evals[i] = num_samples
#                     errors[i] = values[error_type]
#                 ixs = np.argsort(num_evals)
#                 plt.loglog(num_evals[ixs], errors[ixs], "o-",
#                            label=("pce (%s, %s)" % (expansion, sampling_strategy)).replace("_", " "))

        if args.surrogate in ["sg", "both"]:
            for gridType, level, maxGridSize, refinement, isFull in sg_settings:
                key = get_key_sg(gridType, level, maxGridSize, refinement, isFull)
                res = {'t': np.array([]),
                       'mean': np.array([]),
                       'var': np.array([])}
                for t, values in results["sg"][key]["results"].items():
                    res['t'] = np.append(res["t"], t)
                    res['mean'] = np.append(res["mean"], values.values()[-1]["mean_estimated"])
                    res['var'] = np.append(res["var"], values.values()[-1]["var_estimated"])

                ixs = np.argsort(res["t"])
                plt.plot(res['t'][ixs], res[error_type][ixs], "-",
                         linewidth=2,
                         label=("sg (%s, %s)" % (gridType, refinement)).replace("_", " "))


        plt.title(error_type.replace("_", " "))
        plt.legend(loc="lower left")
        plt.show()
