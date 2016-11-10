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
from pysgpp.pysgpp_swig import GridType_Poly, GridType_Linear


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

    def withPolynomialBase(self, deg):
        """
        Defines the polynomial base of the grid
        @param deg: degree of polynomial base as integer
        """
        if deg > 1:
            self.__deg = deg
        else:
            print "Warning: GridDescriptor.withPolynomialBasis - deg < 2 ignored"
        return self

    def withBorder(self, boundaryLevel):
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

    def isClenshawCurtis(self):
        """
        Defines if a clenshaw curtis grid should be generated
        """
        self.__clenshaw_curtis = True
        return self

    def withModifiedBasis(self):
        """
        define a basis with extrapolation towards the boundary
        """
        self.__modified = True
        return self

    def fromGrid(self, grid):
        """
        Indicates that all grid points in grid should also be in the
        new grid
        @param grid:
        """
        self.__grid = grid
        self.__dim = grid.getDimension()
        self.__deg = getDegree(grid)
        if hasBorder(grid):
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

            if self.__grid is not None:
                gridConfig.dim_ = self.__grid.getDimension()
            else:
                gridConfig.dim_ = self.__dim

            if (self.__dim is None or self.level is None) and self.__grid is None:
                raise AttributeError("Not all attributes assigned to create\
                                     grid")
            if self.__boundaryLevel is not None:
                gridConfig.boundaryLevel_ = self.__boundaryLevel

            gridConfig.maxDegree_ = self.__deg

            # identify grid type
            if self.__border is not None:
                if self.__clenshaw_curtis:
                    if self.__deg > 1:
                        gridConfig.type_ = GridType_PolyClenshawCurtisBoundary
                    else:
                        gridConfig.type_ = GridType_LinearClenshawCurtisBoundary
                else:
                    if self.__deg > 1:
                        gridConfig.type_ = GridType_PolyBoundary
                    else:
                        gridConfig.type_ = GridType_LinearBoundary
            else:
                if self.__modified:
                    if self.__clenshaw_curtis:
                        if self.__deg > 1:
                            gridConfig.type_ = GridType_ModPolyClenshawCurtis
                        else:
                            gridConfig.type_ = GridType_ModLinearClenshawCurtis
                    else:
                        if self.__deg > 1:
                            gridConfig.type_ = GridType_ModPoly
                        else:
                            gridConfig.type_ = GridType_ModLinear
                else:
                    if self.__clenshaw_curtis:
                        if self.__deg > 1:
                            gridConfig.type_ = GridType_PolyClenshawCurtis
                        else:
                            gridConfig.type_ = GridType_LinearClenshawCurtis
                    else:
                        if self.__deg > 1:
                            gridConfig.type_ = GridType_Poly
                        else:
                            gridConfig.type_ = GridType_Linear

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
                for i in xrange(copygs.size()):
                    gp = copygs.getPoint(i)
                    # insert grid point
                    if not gs.isContaining(gp):
                        gs.insert(HashGridPoint(gp))
                    if self.__boundaryLevel == 1:
                        insertTruncatedBorder(grid, gp)
                gs.recalcLeafProperty()

        return grid
