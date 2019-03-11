from builtins import range
import os
import pickle as pkl
import numpy as np
from itertools import combinations
from pysgpp.extensions.datadriven.uq.helper import sortPermutations

row_entries = [r"$S_1$ & 0.3138 & 0.3550 & 0.3146 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
               r"$S_2$ & 0.4424 & 0.1846 & 0.4396 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
               r"$S_3$ & 0 & 0.0000 & 0.0000 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
               r"$S_{1, 2}$ & 0 & 0.0000 & 0.0000 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
               r"$S_{1, 3}$ & 0.2436 & 0.4603 & 0.2459 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
               r"$S_{2, 3}$ & 0 & 0.0000 & 0.0000 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f",
               r"$S_{1, 2, 3}$ & 0 & 0.0000 & 0.0000 & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f"]

latexcode = r"""
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
    %s \\
    %s \\
    %s \\
    \hline
    %s
    \bottomrule
  \end{tabularx}
  \caption{blubb}
  \label{tab::ishigami-anova}
\end{table}
"""

unknowns = r"""
    \multicolumn{2}{l}{degree 1d $p$} &
    \multicolumn{1}{l}{$5$} &
    \multicolumn{1}{l}{$9$} &
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
    \multicolumn{1}{l}{$56$} &
    \multicolumn{1}{l}{$220$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\
    \multicolumn{2}{l}{\# model evaluations $N$} &
    \multicolumn{1}{l}{$77$} &
    \multicolumn{1}{l}{$291$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} &
    \multicolumn{1}{l}{$%i$} \\
"""

def get_key_pce(sampling_strategy, degree_1d):
    return (sampling_strategy, degree_1d)

def get_key_sg(gridType, level, maxGridSize, refinement, N):
    return (gridType, level, maxGridSize, refinement, N)

def load_results(path="results"):
    ans = {"pce": {},
           "sg": {}}
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

                if currentStats["surrogate"] == "pce":
                    key = get_key_pce(currentStats["sampling_strategy"],
                                      currentStats["degree_1d"])
                    ans["pce"][key] = currentStats
                elif currentStats["surrogate"] == "sg":
                    key = get_key_sg(currentStats["grid_type"],
                                     currentStats["level"],
                                     currentStats["max_grid_size"],
                                     currentStats["refinement"],
                                     currentStats["results"][-1]["num_model_evaluations"])
                    ans["sg"][key] = currentStats

    return ans


if __name__ == "__main__":
    results = load_results()
    # extract the ones needed for the table
    pce = [('fekete', 5),
           ('fekete', 9),
           ('leja', 5),
           ('leja', 9)]
    sg = [("modPolyClenshawCurtis", 3, 250, None, 31),
          ("modPolyClenshawCurtis", 5, 250, None, 351),
          ("modPolyClenshawCurtis", 2, 65, "var", 73),
          ("modPolyClenshawCurtis", 2, 250, "var", 247)]

    # export values from scenarios to the table
    perms = []
    for k in range(3):
        for perm in combinations([0, 1, 2], r=k + 1):
            perms.append(perm)

    for i, perm in enumerate(sortPermutations(perms)):
        print(i, perm)
        row_entries[i] = row_entries[i] % (np.abs(results["pce"][pce[0]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["pce"][pce[1]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["pce"][pce[2]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["pce"][pce[3]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[0]]["results"][-1]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[1]]["results"][-1]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[2]]["results"][-1]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[3]]["results"][-1]["sobol_indices_estimated"][perm]))

    unknowns = unknowns % (results["pce"][pce[0]]["degree_1d"],
                           results["pce"][pce[1]]["degree_1d"],
                           results["pce"][pce[2]]["degree_1d"],
                           results["pce"][pce[3]]["degree_1d"],
                           results["sg"][sg[0]]["level"],
                           results["sg"][sg[1]]["level"],
                           results["sg"][sg[2]]["level"],
                           results["sg"][sg[3]]["level"],
                           # ---------------------------------------------------------
                           results["pce"][pce[0]]["num_terms"],
                           results["pce"][pce[1]]["num_terms"],
                           results["pce"][pce[2]]["num_terms"],
                           results["pce"][pce[3]]["num_terms"],
                           results["sg"][sg[0]]["results"][-1]["num_model_evaluations"],
                           results["sg"][sg[1]]["results"][-1]["num_model_evaluations"],
                           results["sg"][sg[2]]["results"][-1]["num_model_evaluations"],
                           results["sg"][sg[3]]["results"][-1]["num_model_evaluations"],
                           # ---------------------------------------------------------
                           results["pce"][pce[0]]["num_model_evaluations"],
                           results["pce"][pce[1]]["num_model_evaluations"],
                           results["pce"][pce[2]]["num_model_evaluations"],
                           results["pce"][pce[3]]["num_model_evaluations"],
                           results["sg"][sg[0]]["results"][-1]["num_model_evaluations"],
                           results["sg"][sg[1]]["results"][-1]["num_model_evaluations"],
                           results["sg"][sg[2]]["results"][-1]["num_model_evaluations"],
                           results["sg"][sg[3]]["results"][-1]["num_model_evaluations"]
                           )

    fd = open(os.path.join("tables", "ishigami_results_table.tex"), "w")
    fd.write(latexcode % tuple(row_entries + [unknowns]))
    fd.close()
