import os
import pickle as pkl
import numpy as np

atexcode = r"""
\begin{table}[!ht]
  \fontsize{9pt}{3ex}\selectfont
  \renewcommand{\arraystretch}{1.2}
  \begin{tabularx}{\textwidth}{XXXXXXX}
    \toprule
    index & analytic value &
    \# samples &     \multicolumn{2}{l}{\citet{Sudret2008Global} &     \multicolumn{2}{l}{aPC} &     \multicolumn{2}{l}{SG}\\
    \hline
    & \multicolumn{1}{l}{KL} & \multicolumn{1}{l}{L} &     \multicolumn{1}{l}{KL} & \multicolumn{1}{l}{L} &     \multicolumn{1}{l}{KL} & \multicolumn{1}{l}{L}\\
    \toprule
    50 & 0.2655 & -0.7046 & 2.829 & 1.85 & 0.2157 & -0.7491 \\
    75 & 0.2387 & -0.7314 & 1.837 & 0.858 & 0.1838 & -0.781 \\
    100 & 0.213 & -0.7571 & 1.424 & 0.446 & 0.1533 & -0.8115 \\
    500 & 0.1081 & -0.8655 & 0.4294 & -0.5488 & 0.1157 & -0.849 \\
    1000 & 0.07851 & -0.8951 & 0.2744 & -0.7027 & 0.06948 & -0.8953 \\
    5000 & 0.03964 & -0.9356 & 0.1185 & -0.8598 & 0.02778 & -0.937 \\
    10000 & 0.03001 & -0.9352 & 0.09217 & -0.8847 & 0.02014 & -0.9446 \\
    \bottomrule
  \end{tabularx}
  \caption{}
  \label{tab::ishigami-anova}
\end{table}
"""

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
                    for sampling_strategy in ["leja", "fekete"]:
                        if currentStats["sampling_strategy"] == sampling_strategy:
                            if sampling_strategy not in ans["pce"]:
                                ans["pce"][sampling_strategy] = []
                            ans["pce"][sampling_strategy].append(currentStats)
                elif currentStats["surrogate"] == "sg":
                    for gridType in ["linear", "poly", "modlinear",
                                     "modPolyClenshawCurtis", "modpoly"]:
                        if currentStats["grid_type"] == gridType:
                            if gridType not in ans["sg"]:
                                ans["sg"][gridType] = []
                            ans["sg"][gridType].append(currentStats)

    return ans


if __name__ == "__main__":
    results = load_results()

    # print the pce results
    for sampling_strategy, values in results["pce"].items():
        print sampling_strategy
        for i in np.argsort([values[i]["num_model_evaluations"] for i in xrange(len(values))]):
            print values[i]
