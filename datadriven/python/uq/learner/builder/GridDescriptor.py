from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.learner.formatter.GridFormatter import GridFormatter

from pysgpp import Grid, HashGridIndex

import os
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import (insertTruncatedBorder,
                                           hasBorder,
                                           getDegree)


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
        self.__deg = None
        self.level = None
        self.__file = None
        self.__border = None
        self.__grid = None
        self.__full = None

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
        self.__deg = deg
        return self

    def withBorder(self, border):
        """
        Defines the border type of the grid
        @param border: border type as defined in bin.learner.Types.BorderTypes
        """
        self.__border = border
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
        self.__dim = grid.getStorage().dim()
        self.__deg = getDegree(grid)
        self.__border = hasBorder(grid)
        if self.__border:
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
            if self.__grid is not None:
                self.__dim = self.__grid.getStorage().dim()

            if (self.__dim is None or self.level is None) and self.__grid is None:
                raise AttributeError("Not all attributes assigned to create\
                                     grid")
            if self.__border is not None:
                if self.__border == BorderTypes.TRAPEZOIDBOUNDARY:
                    if self.__deg > 1:
                        grid = Grid.createPolyTruncatedBoundaryGrid(self.__dim, self.__deg)
                    else:
                        grid = Grid.createLinearTruncatedBoundaryGrid(self.__dim)
                elif self.__border == BorderTypes.COMPLETEBOUNDARY:
                    if self.__deg > 1:
                        raise NotImplementedError()
                    else:
                        grid = Grid.createLinearBoundaryGrid(self.__dim)
                else:
                    if self.__deg > 1:
                        grid = Grid.createModPolyGrid(self.__dim, self.__deg)
                    else:
                        grid = Grid.createModLinearGrid(self.__dim)
            else:
                # no border points
                if self.__deg > 1:
                    grid = Grid.createPolyGrid(self.__dim, self.__deg)
                else:
                    grid = Grid.createLinearGrid(self.__dim)

            # generate the grid
            if self.level is not None:
                generator = grid.createGridGenerator()
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
                    gp = copygs.get(i)
                    # insert grid point
                    if not gs.has_key(gp):
                        gs.insert(HashGridIndex(gp))
                    if self.__border == BorderTypes.TRAPEZOIDBOUNDARY:
                        insertTruncatedBorder(grid, gp)
                gs.recalcLeafProperty()

        return grid
