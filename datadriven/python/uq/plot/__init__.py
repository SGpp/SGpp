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
                    plotHistogram1d,
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
                    plotTimedependentDensity2dWithRawData,
                    plotTimedependentDensity2d,
                    plotDensity2d,
                    plotSG2d,
                    plotSamples2d,
                    plotGrid2d,
                    plotFunction2d,
                    plotScatter2d)

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
    from colors import initialize_plotting_style

    initialize_plotting_style()

    __all__.append('plt')
    __all__.append('mpl')
except:
    print "Error initializing plotting style"
    pass
