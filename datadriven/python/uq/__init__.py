"""
UQ package for SGpp
===============================
"""

__version__ = "1.0"

__all__ = ["analysis",
           "dists",
           "estimators",
           "learner",
           "operations",
           "parameters",
           "plot",
           "quadrature",
           "refinement",
           "sampler",
           "transformation",
           "uq_setting",
           "manager",
           "jsonLib",
           "tools",
           "toolsKbhitCountdown"]

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

try:
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('font', **{'family':'sans-serif',
                  'sans-serif':['Helvetica'],
                  'size': 20})
    rc('text', usetex=True)
    # must be included in file set options do not seem to be exported
    # rc('xtick', labelsize=20)
    # rc('ytick', labelsize=20)


    import matplotlib as mpl
    # Use true LaTeX and bigger font
    mpl.rc('text', usetex=True)
    # Include packages `amssymb` and `amsmath` in LaTeX preamble
    # as they include extended math support (symbols, envisonments etc.)
    mpl.rcParams['text.latex.preamble'] = [r"\usepackage[utf8]{inputenc}",
                                           r"\usepackage{tikz}",
                                           r"\usepackage{pgfplots}",
                                           r"\usepackage{amssymb}",
                                           r"\usepackage{amsmath}", ]
    mpl.rcParams['axes.labelsize'] = 20
    __all__.append('plt')
    __all__.append('rc')
    __all__.append('mpl')
except:
    pass
