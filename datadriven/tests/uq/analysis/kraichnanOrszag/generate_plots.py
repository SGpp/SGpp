from builtins import str
import os
import pickle as pkl
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.helper import sortPermutations
from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend, \
    load_color, load_marker, savefig


def get_key_pce(expansion, sampling_strategy, N):
    return (expansion, sampling_strategy, N)

def get_key_sg(gridType, maxGridSize, refinement, isFull):
    return (gridType, maxGridSize, refinement, isFull)


def get_key_mc(sampling_strategy, N):
    return (sampling_strategy, N)


def parse_monte_carlo_results(results):
    key = get_key_mc(sampling_strategy, numSamples)
    time_steps = results["mc"][key]["time_steps"]
    res = {}
    for t in np.sort(time_steps):
        res[t] = {}
        for key_value in ["mean", "var"]:
            key_error = "%s_error" % key_value
            values = results["mc"][key]["results"]["%s_estimated" % key_value][t]

            res[t][key_value] = values["value"]
            res[t][key_error] = values["confidence_interval"]

    return res


def parse_sg_results(gridType,
                     maxGridSize,
                     refinement,
                     isFull,
                     results):
    key = get_key_sg(gridType, maxGridSize, refinement, isFull)
    time_steps = results["sg"][key]["time_steps"]
    for t in np.sort(time_steps):
        res[t] = {}
        for key_value in ["mean", "var"]:
            values = results["sg"][key]["results"][t]
            maxLevel = max(values.keys())
            res[t][key_value] = values[maxLevel]["%s_estimated_per_iteration" % key_value][0]
            res[t]["sobol_indices"] = values[maxLevel]["sobol_indices"]
    return res


def load_results(inputspace, setting, qoi, settings, path="results"):
    ans = {"pce": {},
           "sg": {},
           "mc": {}}
    for root, dirs, files in os.walk(os.path.join(path, inputspace)):
        for filename in files:
            if "pkl" in filename and \
                    "qoi%s" % qoi in filename and \
                    "s%i" % setting in filename:
                load_file = True
                if "sg" in filename:
                    load_file = False
                    for _, maxGridSize, _, _ in settings[setting]["sg"][qoi]:
                        if str(maxGridSize) in filename:
                            load_file = True
                            break

                if load_file:
                    print("-" * 80)
                    print("load %s -> " % filename)
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
                        elif "sg" in currentStats["surrogate"]:
                            key = get_key_sg(currentStats["grid_type"],
                                             currentStats["max_grid_size"],
                                             currentStats["refinement"],
                                             currentStats["is_full"])

                            ans["sg"][key] = currentStats
                        else:
                            key = get_key_mc(currentStats["sampling_strategy"],
                                             currentStats["num_model_evaluations"])
                            ans["mc"][key] = currentStats

                        print("  %s" % (key,))
                else:
                    print("-" * 80)
                    print("available %s -> " % filename)
    return ans


def plotMCResults(ts, values, err, color, marker, label):
    plt.plot(ts, values, label=label, marker=marker, color=color)
    plt.fill_between(ts, values, err[:, 0],
                     facecolor=color, alpha=0.2)
    plt.fill_between(ts, err[:, 1], values,
                     facecolor=color, alpha=0.2)


settings = {1: {'mc': [('latin_hypercube', 2000)],
                'sg': [
                       (8, 1025, None, False),
                       (8, 617, "simple", False)
                       ],
                'pce': [("full_tensor", 'gauss', 4000),
                        ('total_degree', 'gauss_leja', 4000)]},
            2: {'mc': [('latin_hypercube', 2000)],
                'sg': [
                       (8, 1281, None, False),
                       (8, 1023, 'simple', False)
                       ],
                'pce': [("full_tensor", 'gauss', 4000),
                        ('total_degree', 'gauss_leja', 4000)]},
            3: {'mc': [('latin_hypercube', 100000)],
                'sg': {"y1": [
                              (8, 1505, None, False),
                              (3, 783, 'simple', False),
                              (8, 643, 'simple', False),
                              (8, 689, 'l2', False),
                              (3, 757, 'l2', False)
                              ],
                       "y3":[
                             (8, 1505, None, False),
                             (3, 793, 'simple', False),
                             (8, 821, 'simple', False),
                             (8, 745, 'l2', False),
                             (3, 701, 'l2', False)
                             ]},
                'pce': [("full_tensor", 'gauss', 4000),
                        ('total_degree', 'gauss_leja', 4000)]}}

if __name__ == "__main__":
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="uniform", type=str, help="define which probabilistic model should be used")
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--setting', default=3, type=int, help='parameter settign for test problem')
    parser.add_argument('--qoi', default="y1", type=str, help="define the quantity of interest")
    args = parser.parse_args()

    results = load_results(args.model, args.setting, args.qoi, settings)

    error_type = "var"
    # extract the ones needed for the table
    mc_settings = settings[args.setting]["mc"]
    pce_settings = settings[args.setting]["pce"]
    sg_settings = settings[args.setting]["sg"][args.qoi]

    # plot mc results to compare
    fig = plt.figure()
    for sampling_strategy, numSamples in mc_settings:
        res = parse_monte_carlo_results(results)
        time_steps = np.array(list(res.keys()))
        ixs = np.argsort(time_steps)
        time_steps = time_steps[ixs]
        ixs = np.where(time_steps <= 6)[0]
        time_steps = time_steps[ixs]
        values = np.ndarray(time_steps.shape)
        err = np.ndarray((time_steps.size, 2))
        for i, t in enumerate(time_steps[ixs]):
            values[i] = res[t][error_type]
            err[i, :] = res[t]["%s_error" % error_type]

        plotMCResults(time_steps[ixs],
                      values,
                      err,
                      color=load_color(0),
                      marker=load_marker(0),
                      label=r"MC (M=$10^5$)")

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
        for k, (gridType, maxGridSize, refinement, isFull) in enumerate(sg_settings):
            res = parse_sg_results(gridType,
                                   maxGridSize,
                                   refinement,
                                   isFull,
                                   results)
            time_steps = np.array(list(res.keys()))
            ixs = np.argsort(time_steps)
            time_steps = time_steps[ixs]
            ixs = np.where(time_steps <= 6)[0]
            time_steps = time_steps[ixs]
            values = np.ndarray(time_steps.shape)
            err = np.ndarray((time_steps.size, 2))
            for i, t in enumerate(time_steps):
                values[i] = res[t][error_type]
                err[i, :] = values[i]

            if refinement is None:
                label = r"SG ($N=%i$)" % (maxGridSize,)
            else:
                label = r"aSG ($N=%i$)" % (maxGridSize,)

            plotMCResults(time_steps,
                          values,
                          err,
                          color=load_color(k + 1),
                          marker=load_marker(k + 1),
                          label=label)

    plt.ylabel("$\mathbb{V}$ of $y_%s$" % args.qoi[-1])
    lgd = insert_legend(fig, loc="bottom", ncol=2)
    savefig(fig, "plots/kraichnan_orszag_%s_s%i_%s_%s" % (args.model,
                                                          args.setting,
                                                          args.qoi,
                                                          error_type))
    plt.close(fig)
