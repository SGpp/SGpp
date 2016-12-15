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
    \multicolumn{4}{c}{with modfied polynomial basis} \\
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
    \multicolumn{2}{c}{degree 1d $p$} &
    \multicolumn{1}{c}{$5$} &
    \multicolumn{1}{c}{$9$} &
    \multicolumn{1}{c}{$5$} &
    \multicolumn{1}{c}{$9$} &
    \multicolumn{1}{c}{$5$} &
    \multicolumn{1}{c}{$9$} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} \\
    \multicolumn{2}{c}{grid level $\level$} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{-} &
    \multicolumn{1}{c}{$3$} &
    \multicolumn{1}{c}{$5$} &
    \multicolumn{1}{c}{$2$} &
    \multicolumn{1}{c}{$2$} \\
    \multicolumn{2}{c}{\# unknown coefficients} &
    \multicolumn{1}{c}{$56$} &
    \multicolumn{1}{c}{$220$} &
    \multicolumn{1}{c}{$56$} &
    \multicolumn{1}{c}{$220$} &
    \multicolumn{1}{c}{$56$} &
    \multicolumn{1}{c}{$220$} &
    \multicolumn{1}{c}{$31$} &
    \multicolumn{1}{c}{$351$} &
    \multicolumn{1}{c}{$73$} &
    \multicolumn{1}{c}{$247$} \\
    \multicolumn{2}{c}{\# model evaluations $N$} &
    \multicolumn{1}{c}{$77$} &
    \multicolumn{1}{c}{$291$} &
    \multicolumn{1}{c}{$67$} &
    \multicolumn{1}{c}{$264$} &
    \multicolumn{1}{c}{$67$} &
    \multicolumn{1}{c}{$264$} &
    \multicolumn{1}{c}{$31$} &
    \multicolumn{1}{c}{$351$} &
    \multicolumn{1}{c}{$73$} &
    \multicolumn{1}{c}{$247$} \\
    \bottomrule
  \end{tabularx}
  \caption{blubb}
  \label{tab::ishigami-anova}
\end{table}

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
                print "=" * 80
                print "loading '%s'" % path
                print "=" * 80
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
    for k in xrange(3):
        for perm in combinations([0, 1, 2], r=k + 1):
            perms.append(perm)
            
    for i, perm in enumerate(sortPermutations(perms)):
        print i, perm
        row_entries[i] = row_entries[i] % (np.abs(results["pce"][pce[0]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["pce"][pce[1]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["pce"][pce[2]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["pce"][pce[3]]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[0]]["results"][-1]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[1]]["results"][-1]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[2]]["results"][-1]["sobol_indices_estimated"][perm]),
                                           np.abs(results["sg"][sg[3]]["results"][-1]["sobol_indices_estimated"][perm]))

    fd = open(os.path.join("tables", "ishigami_results_table.tex"), "w")
    fd.write(latexcode % tuple(row_entries))
    fd.close()
