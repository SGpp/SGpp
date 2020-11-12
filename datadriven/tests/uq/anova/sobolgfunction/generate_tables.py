# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import os
import pickle as pkl
import numpy as np
from itertools import combinations
from pysgpp.extensions.datadriven.uq.helper import sortPermutations

latexcode_full_model = r"""
\begin{table}[!ht]
  \fontsize{8pt}{3ex}\selectfont
  \renewcommand{\arraystretch}{1.2}
  \begin{tabularx}{\textwidth}{XXXXXX}
    \toprule
    index &
    analytic &
    \multicolumn{2}{c}{generalized PC expansion} &
    \multicolumn{2}{c}{sparse grid on Clenshaw-Curtis} \\
    & value & & &
    \multicolumn{2}{c}{with modified polynomial basis} \\
    \hline
    & &
    \multicolumn{1}{l}{Fekete} &
    \multicolumn{1}{l}{Leja} &
    \multicolumn{1}{l}{regular} &
    \multicolumn{1}{l}{adaptive} \\
    \toprule
    %s \\
    %s \\
    %s \\
    %s \\
    %s \\
    %s \\
    %s \\
    %s \\
    \hline
    %s
    \bottomrule
  \end{tabularx}
  \caption{blubb}
  \label{tab::sobolgfunction-full-model-anova}
\end{table}
"""

latexcode_reduced_model = r"""
\begin{table}[!ht]
  \fontsize{8pt}{3ex}\selectfont
  \renewcommand{\arraystretch}{1.2}
  \begin{tabularx}{\textwidth}{XXXXXXXXXXXX}
    \toprule
    index &
    analytic &
    \multicolumn{6}{c}{generalized PC expansion} &
    \multicolumn{4}{c}{sparse grid on Clenshaw-Curtis} \\
    &
    value &
    & & & & & &
    \multicolumn{4}{c}{with modified polynomial basis} \\
    \hline
    & &
    \multicolumn{2}{c}{\citet{Sudret08Global}} &
    \multicolumn{2}{c}{Fekete} &
    \multicolumn{2}{c}{Leja} &
    \multicolumn{2}{c}{regular} &
    \multicolumn{2}{c}{adaptive} \\
    \toprule
    %s \\
    %s \\
    %s \\
    %s \\
    \hline
    %s
    \bottomrule
  \end{tabularx}
  \caption{blubb}
  \label{tab::sobolgfunction-reduced-model-anova}
\end{table}
"""


row_entries = {"full": [r"$ST_1$ & 0.6342 & %.4f & %.4f & %.4f & %.4f",
                        r"$ST_2$ & 0.2683 & %.4f & %.4f & %.4f & %.4f",
                        r"$ST_3$ & 0.0671 & %.4f & %.4f & %.4f & %.4f",
                        r"$ST_4$ & 0.0200 & %.4f & %.4f & %.4f & %.4f",
                        r"$ST_5$ & 0.0055 & %.4f & %.4f & %.4f & %.4f",
                        r"$ST_6$ & 0.0009 & %.4f & %.4f & %.4f & %.4f",
                        r"$ST_7$ & 0.0002 & %.4f & %.4f & %.4f & %.4f",
                        r"$ST_8$ & 0.0000 & %.4f & %.4f & %.4f & %.4f"],
               "reduced": [r"$ST_1$ & 0.6342 & 0.6644 & 0.6350 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
                           r"$ST_2$ & 0.2945 & 0.2611 & 0.3057 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
                           r"$ST_3$ & 0.0756 & 0.0581 & 0.0756 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
                           r"$ST_4$ & 0.0227 & 0.0164 & 0.0220 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f"]}


unknowns = {"full": r"""    \multicolumn{2}{l}{degree 1d $p$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} \\
    \multicolumn{2}{l}{grid level $\level$} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\
    \multicolumn{2}{l}{\# unknown coefficients} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\
    \multicolumn{2}{l}{\# model evaluations $N$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\""",
                        "reduced": r"""    \multicolumn{2}{l}{degree 1d $p$} &
    \multicolumn{1}{l}{$3$} &
    \multicolumn{1}{l}{$5$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} \\
    \multicolumn{2}{l}{grid level $\level$} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{-} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\
    \multicolumn{2}{l}{\# unknown coefficients} &
    \multicolumn{1}{l}{$35$} &
    \multicolumn{1}{l}{$126$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\
    \multicolumn{2}{l}{\# model evaluations $N$} &
    \multicolumn{1}{l}{$73$} &
    \multicolumn{1}{l}{$233$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\"""}

latexcode = {"full": latexcode_full_model,
             "reduced": latexcode_reduced_model}

def get_key_pce(sampling_strategy, degree_1d, num_samples):
    return (sampling_strategy, degree_1d, num_samples)

def get_key_sg(gridType, level, maxGridSize, refinement, N):
    return (gridType, level, maxGridSize, refinement, N)

def load_results(path="results"):
    ans = {"full": {"sg": {},
                    "pce": {}},
           "reduced": {"pce": {},
                       "sg": {}}}
    for root, dirs, files in os.walk(path):
        for filename in files:
            if "pkl" in filename:
                path = os.path.join(root, filename)
                print("=" * 80)
                print("loading '%s'" % path)
                print("=" * 80)
                fd = open(path, "r")
                currentStats = pkl.load(fd)
                fd.close()

                if currentStats["num_dims"] == 4:
                    model = "reduced"
                else:
                    model = "full"

                if currentStats["surrogate"] == "pce":
                    key = get_key_pce(currentStats["sampling_strategy"],
                                      currentStats["degree_1d"],
                                      currentStats["num_model_evaluations"])
                    ans[model]["pce"][key] = currentStats
                elif currentStats["surrogate"] == "sg":
                    key = get_key_sg(currentStats["grid_type"],
                                     currentStats["level"],
                                     currentStats["max_grid_size"],
                                     currentStats["refinement"],
                                     currentStats["results"][-1]["num_model_evaluations"])
                    ans[model]["sg"][key] = currentStats

    return ans


if __name__ == "__main__":
    results = load_results()
    # extract the ones needed for the table
    pce = {'full': [('fekete', 2, 72),
                    ('leja', 2, 72)],
           'reduced': [('fekete', 3, 56),
                       ('fekete', 5, 201),
                       ('leja', 3, 56),
                       ('leja', 5, 201)]}
    sg = {'full': [("modPolyClenshawCurtis", 3, 140, None, 161),
                   ("modPolyClenshawCurtis", 2, 70, "var", 75)],
          'reduced': [("modPolyClenshawCurtis", 3, 200, None, 49),
                      ("modPolyClenshawCurtis", 4, 200, None, 209),
                      ("modPolyClenshawCurtis", 2, 100, "var", 103),
                      ("modPolyClenshawCurtis", 2, 200, "var", 201)]}

    # export values from scenarios to the table
    model = "full"
    for i, perm in enumerate(range(8)):
        perm = (perm,)
        row_entries[model][i] = row_entries[model][i] % (np.abs(results[model]["pce"][pce[model][0]]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["pce"][pce[model][1]]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["sg"][sg[model][0]]["results"][-1]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["sg"][sg[model][1]]["results"][-1]["total_effects_estimated"][perm]))

    unknowns[model] = unknowns[model] % (results[model]["pce"][pce[model][0]]["degree_1d"],
                                         results[model]["pce"][pce[model][1]]["degree_1d"],
                                         results[model]["sg"][sg[model][0]]["level"],
                                         results[model]["sg"][sg[model][1]]["level"],
                                         # ---------------------------------------------------------
                                         results[model]["pce"][pce[model][0]]["num_terms"],
                                         results[model]["pce"][pce[model][1]]["num_terms"],
                                         results[model]["sg"][sg[model][0]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][1]]["results"][-1]["num_model_evaluations"],
                                         # ---------------------------------------------------------
                                         results[model]["pce"][pce[model][0]]["num_model_evaluations"],
                                         results[model]["pce"][pce[model][1]]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][0]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][1]]["results"][-1]["num_model_evaluations"]
                                         )

    fd = open(os.path.join("tables", "sobolgfunction_%s_results_table.tex" % model), "w")
    fd.write(latexcode[model] % tuple(row_entries[model] + [unknowns[model]]))
    fd.close()

    model = "reduced"
    for i, perm in enumerate(range(4)):
        perm = (perm,)
        row_entries[model][i] = row_entries[model][i] % (np.abs(results[model]["pce"][pce[model][0]]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["pce"][pce[model][1]]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["pce"][pce[model][2]]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["pce"][pce[model][3]]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["sg"][sg[model][0]]["results"][-1]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["sg"][sg[model][1]]["results"][-1]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["sg"][sg[model][2]]["results"][-1]["total_effects_estimated"][perm]),
                                                         np.abs(results[model]["sg"][sg[model][3]]["results"][-1]["total_effects_estimated"][perm]))

    unknowns[model] = unknowns[model] % (results[model]["pce"][pce[model][0]]["degree_1d"],
                                         results[model]["pce"][pce[model][1]]["degree_1d"],
                                         results[model]["pce"][pce[model][2]]["degree_1d"],
                                         results[model]["pce"][pce[model][3]]["degree_1d"],
                                         results[model]["sg"][sg[model][0]]["level"],
                                         results[model]["sg"][sg[model][1]]["level"],
                                         results[model]["sg"][sg[model][2]]["level"],
                                         results[model]["sg"][sg[model][3]]["level"],
                                         # ---------------------------------------------------------
                                         results[model]["pce"][pce[model][0]]["num_terms"],
                                         results[model]["pce"][pce[model][1]]["num_terms"],
                                         results[model]["pce"][pce[model][2]]["num_terms"],
                                         results[model]["pce"][pce[model][3]]["num_terms"],
                                         results[model]["sg"][sg[model][0]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][1]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][2]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][3]]["results"][-1]["num_model_evaluations"],
                                         # ---------------------------------------------------------
                                         results[model]["pce"][pce[model][0]]["num_model_evaluations"],
                                         results[model]["pce"][pce[model][1]]["num_model_evaluations"],
                                         results[model]["pce"][pce[model][2]]["num_model_evaluations"],
                                         results[model]["pce"][pce[model][3]]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][0]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][1]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][2]]["results"][-1]["num_model_evaluations"],
                                         results[model]["sg"][sg[model][3]]["results"][-1]["num_model_evaluations"]
                                         )

    fd = open(os.path.join("tables", "sobolgfunction_%s_results_table.tex" % model), "w")
    fd.write(latexcode[model] % tuple(row_entries[model] + [unknowns[model]]))
    fd.close()
