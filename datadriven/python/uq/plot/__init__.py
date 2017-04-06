"""
plot package
==========================================

"""
__version__ = "0.1"

__all__ = ["plotSGDE3d", "plotSGDE2d",
           "plotGrid"]

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"


from plot1d import (plotFunction1d,
                    plotDensity1d,
                    plotCDF,
                    plotCDF1d,
                    plotPDF,
                    plotGrid1d,
                    plotSG1d,
                    plotSurplusLevelWise,
                    plotNodal1d,
                    plotSGNodal1d,
                    plotSGDE1d,
                    plotSobolIndices)

from plot2d import (plotSGDE2d,
                    plotDensity2d,
                    plotSG2d,
                    plotSamples2d,
                    plotGrid2d,
                    plotFunction2d)

from plot3d import (plotFunction3d,
                    plotDensity3d,
                    plotKDE3d,
                    plotSG3d,
                    plotGrid3d,
                    plotNodal3d,
                    plotSGNodal3d,
                    plotError3d,
                    insert_labels_3d)

from plotGrid import plotGrid

from scatterplot import scatterplot_matrix

from colors import intToRGB, rgbTpInt, load_color, load_marker, load_font, \
    load_font_properties, insert_legend

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from colors import load_font, load_default_color_map
    from mpltools import layout

    pysgpp_uq_font = load_font()

    plt.style.use('seaborn-paper')

    # Include packages `amssymb` and `amsmath` in LaTeX preamble
    # as they include extended math support (symbols, envisonments etc.)
    pgf_with_custom_preamble = {"font.family": pysgpp_uq_font["family"],  # use serif/main font for text elements
                                "text.usetex": True,  # use inline math for ticks
                                "text.latex.preamble": [
                                                 r'\usepackage{amsmath}',
                                                 r'\usepackage[scientific-notation=true]{siunitx}',
                                                 r"\usepackage[utf8x]{inputenc}",
                                                 r"\usepackage[T1]{fontenc}",
                                                 r"\usepackage{tikz}",
                                                 r"\usepackage{pgfplots}",
                                                 r"\usepackage{amssymb}",
                                                 r"\usepackage{amsmath}"
                                                 ],
                                'axes.labelsize': pysgpp_uq_font["size"],
                                'font.size': pysgpp_uq_font["size"],
                                'legend.fontsize': pysgpp_uq_font["size"],
                                'xtick.labelsize': pysgpp_uq_font["size"],
                                'ytick.labelsize': pysgpp_uq_font["size"],
                                'axes.unicode_minus': True,
                                'figure.figsize': (5, 4.5),
                                'image.cmap': load_default_color_map(dtype="string")
                                }
    mpl.rcParams.update(pgf_with_custom_preamble)

    __all__.append('plt')
    __all__.append('rc')
    __all__.append('mpl')
except:
    print "Error"
    pass
