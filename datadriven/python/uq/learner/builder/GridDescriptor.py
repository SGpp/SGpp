from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.learner.formatter.GridFormatter import GridFormatter

from pysgpp import Grid, HashGridPoint

import os
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import (insertTruncatedBorder,
                                           hasBorder,
                                           getDegree)
from pysgpp import RegularGridConfiguration, GridType_PolyBoundary, \
    GridType_PolyClenshawCurtis, GridType_PolyClenshawCurtisBoundary, \
    GridType_LinearClenshawCurtisBoundary, GridType_LinearClenshawCurtis, \
    GridType_ModPolyClenshawCurtis, GridType_ModLinearClenshawCurtis, \
    GridType_ModLinear, GridType_ModPoly, GridType_LinearBoundary
from pysgpp.pysgpp_swig import GridType_Poly, GridType_Linear, GridType_Bspline, \
    GridType_BsplineBoundary, GridType_BsplineClenshawCurtis, \
    GridType_ModBsplineClenshawCurtis, GridType_ModBspline


class GridDescriptor(object):
    """
    Grid Descriptor helps to implement fluid interface patter on python
    it encapsulates functionality concerning creation of the grid
    """

    def __init__(self):
        """
        Constructor
        """
        self.__dim = None
        self.__deg = 1
        self.level = None
        self.__file = None
        self.__boundaryLevel = None
        self.__grid = None
        self.__gridType = GridType_Linear
        self.__full = None
        self.__clenshaw_curtis = False
        self.__modified = False

    def withDimension(self, dim):
        """
        Defines the dimensionality of the grid
        @param dim: dimensionality as integer
        """
        self.__dim = dim
        return self

    def withLevel(self, level):
        """
        Defines the level of the grid
        @param level: level as integer
        """
        self.level = level
        return self

    def hasType(self, gridType):
        """
        Defines the grid type
        @param gridType: type of the grid
        """
        self.__gridType = gridType
        return self

    def withDegree(self, deg):
        """
        Defines the polynomial degree of the basis
        @param deg: degree of polynomial base as integer
        """
        if deg > 1:
            self.__deg = deg
        else:
            print "Warning: GridDescriptor.withDegree - deg < 2 ignored"
        return self

    def withBoundaryLevel(self, boundaryLevel):
        """
        Defines the border type of the grid
        @param boundaryLevel: level of the boundary
        """
        self.__boundaryLevel = boundaryLevel
        return self

    def isFull(self):
        """
        Defines if a full grid should be generated
        """
        self.__full = True
        return self

    def fromGrid(self, grid):
        """
        Indicates that all grid points in grid should also be in the
        new grid
        @param grid:
        """
        self.__grid = grid
        self.__dim = grid.getStorage().getDimension()
        self.__deg = getDegree(grid)
        self.__gridType = grid.getType()
        if hasBorder(grid.getType()):
            self.__boundaryLevel = 1
            self.level = 0

        return self

    def fromFile(self, filename):
        """
        Indicates that grid should be restored from file
        @param filename: string name of file the grid should be restored from
        """
        self.__file = filename
        return self

    def createGrid(self):
        """
        Creates the specified grid
        """
        grid = None
        if self.__file is not None and os.path.exists(self.__file):
            gridFormatter = GridFormatter()
            grid = gridFormatter.deserializeFromFile(self.__file)
        else:
            gridConfig = RegularGridConfiguration()
            gridConfig.dim_ = self.__dim

            if (self.__dim is None or self.level is None) and self.__grid is None:
                raise AttributeError("Not all attributes assigned to create\
                                     grid")
            if self.__boundaryLevel is not None:
                gridConfig.boundaryLevel_ = self.__boundaryLevel

            gridConfig.maxDegree_ = self.__deg

            if self.__gridType not in [GridType_Linear,
                                       GridType_LinearBoundary,
                                       GridType_ModLinear,
                                       GridType_LinearClenshawCurtis,
                                       GridType_LinearClenshawCurtisBoundary,
                                       GridType_ModLinearClenshawCurtis,
                                       GridType_Poly,
                                       GridType_PolyBoundary,
                                       GridType_ModPoly,
                                       GridType_PolyClenshawCurtis,
                                       GridType_PolyClenshawCurtisBoundary,
                                       GridType_ModPolyClenshawCurtis,
                                       GridType_Bspline,
                                       GridType_ModBspline,
                                       GridType_BsplineBoundary,
                                       GridType_BsplineClenshawCurtis,
                                       GridType_ModBsplineClenshawCurtis]:
                print "Warning: grid type not fully supported"

            gridConfig.type_ = self.__gridType

            # generate the grid
            grid = Grid.createGrid(gridConfig)

            if self.level is not None:
                generator = grid.getGenerator()
                if not self.__full:
                    generator.regular(self.level)
                else:
                    generator.full(self.level)

            # if there is a grid specified, add all the missing points
            if self.__grid is not None:
                gs = grid.getStorage()
                copygs = self.__grid.getStorage()

                # insert grid points
                for i in xrange(copygs.getSize()):
                    gp = copygs.getPoint(i)
                    # insert grid point
                    if not gs.isContaining(gp):
                        gs.insert(HashGridPoint(gp))
                    if self.__boundaryLevel == 1:
                        insertTruncatedBorder(grid, gp)
                gs.recalcLeafProperty()

        return grid
