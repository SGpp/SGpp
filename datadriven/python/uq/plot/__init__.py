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
                    plotNodal3d,
                    plotSGNodal3d)

from plotGrid import plotGrid

from scatterplot import scatterplot_matrix

from colors import intToRGB, rgbTpInt, loadColorSequence

