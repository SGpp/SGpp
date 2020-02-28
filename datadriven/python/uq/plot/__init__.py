# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
plot package
==========================================

"""
__version__ = "0.1"

__all__ = ["plotSGDE3d", "plotSGDE2d",
           "plotGrid"]

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"


from pysgpp.extensions.datadriven.uq.plot.plot1d import (plotFunction1d,
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

from pysgpp.extensions.datadriven.uq.plot.plot2d import (plotSGDE2d,
                    plotTimedependentDensity2dWithRawData,
                    plotTimedependentDensity2d,
                    plotDensity2d,
                    plotSG2d,
                    plotSamples2d,
                    plotGrid2d,
                    plotFunction2d,
                    plotScatter2d)

from pysgpp.extensions.datadriven.uq.plot.plot3d import (plotFunction3d,
                    plotDensity3d,
                    plotKDE3d,
                    plotSG3d,
                    plotGrid3d,
                    plotNodal3d,
                    plotSGNodal3d,
                    plotError3d,
                    insert_labels_3d)

from pysgpp.extensions.datadriven.uq.plot.plotGrid import plotGrid

from pysgpp.extensions.datadriven.uq.plot.scatterplot import scatterplot_matrix

from pysgpp.extensions.datadriven.uq.plot.colors import intToRGB, rgbTpInt, load_color, load_marker, load_font, \
    load_font_properties, insert_legend

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from pysgpp.extensions.datadriven.uq.plot.colors import initialize_plotting_style

    initialize_plotting_style()

    __all__.append('plt')
    __all__.append('mpl')
except:
    print( "Error initializing plotting style" )
    pass
