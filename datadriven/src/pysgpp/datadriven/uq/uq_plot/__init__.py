"""
plot package
==========================================

"""

__version__ = "0.1"

__all__ = ["plotSGDE3d", "plotSGDE2d",
           "plotGrid"]

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from plot1d import (plotDensity1d,
                    plotCDF,
                    plotPDF,
                    plotSG1d,
                    plotSurplusLevelWise,
                    plotNodal1d,
                    plotSGDE1d)

from plot2d import (plotSGDE2d,
                    plotDensity2d,
                    plotSG2d,
                    plotSamples2d,
                    plotGrid2d,
                    plotFunction2d)

from plot3d import (plot3d,
                    plotDensity3d,
                    plotKDE3d,
                    plotSG3d,
                    plotNodal3d,
                    plotSGNodal3d)

from plotGrid import plotGrid

from scatterplot import scatterplot_matrix
