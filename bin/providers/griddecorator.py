# This file is part of sgpp, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.
#

from bin.pysgpp import *
from bin.providers.evalproviders import *
from bin.providers.refineproviders import *

##
#
class GridDecorator(object):
    """Grid status type. Contains everything needed for datamining"""

    def __init__(self, grid, mode, basetype, zeh):
        self.grid = grid
        self.mode = mode
        self.basetype = basetype
        self.zeh = zeh
        
        self.__refine = None
        self.__eval = None

    def getRefine(self):
        return self.__refine


#    def setRefine(self, value):
#        self.__refine = value
    ##
    # @param type String 
    def setRefine(self, type):
        self.__refine = self.refine_providers[type](self)

    def getEval(self):
        return self.__eval

#    def setEval(self, value):
#        self.__eval = value
    
    ##
    # @param type String 
    def setEval(self, type):
        self.__eval = self.eval_providers[type](self)

    def __getstate__(self):
        odict = self.__dict__.copy()
        del odict['grid']
        
        del odict['_GridDecorator__refine']
        del odict['_GridDecorator__eval']
        
        odict['grid'] = self.grid.serialize()
        
        return odict
    
    def __setstate__(self, dict):
        #restore grid
        self.__dict__.update(dict)
        self.grid = Grid.unserialize(dict['grid'])
    
    def getSize(self):
        return self.grid.getStorage().size()   
        
    ## Constructs GridDecorator object
    # Deserializes a GridDecorator object from options.grid file or creates a new GridDecorator object
    # @return GridDecorator object
    @classmethod
    def create(cls, options, dim = None, mode = None):
        gridDecorator = None
        if options.grid:
            import pickle
            fin = open(options.grid, "rb")
            gridDecorator = pickle.load(fin)
            fin.close()
        else:
            gridDecorator = GridDecorator(cls.constructGrid(options, dim), mode, options.basetype, options.zeh)
        
        return gridDecorator
    
    ## Constructs a new grid.
    # If options.grid is set, then read in a stored grid. If not, construct a new
    # grid dependent on the dimension dim, on options.level and options.polynom.
    # Sets the use of boundary functions according to options.border.
    # @param dim the grid dimension
    # @return a grid
    @classmethod
    def constructGrid(cls, options, dim):
        factories = {"linear" : lambda : Grid.createLinearGrid(dim),
                 "modlinear" : lambda : Grid.createModLinearGrid(dim),
                 "poly" : lambda : Grid.createPolyGrid(dim, options.polynom),
                 "modpoly" : lambda : Grid.createModPolyGrid(dim, options.polynom),
                 }
        grid = factories[options.basetype]()
        gen = grid.createGridGenerator()
        gen.regular(options.level)
        
        return grid
    
    ## List of available refine providers
    # See ClassesEvalProvider for an example        
    eval_providers = {
            "classes" : ClassesEvalProvider,
            "regression" : RegressionEvalProvider,
            }
    
    ##list of refine providers
    # see SuprlusRefineProvider for an example
    refine_providers = {
            "surplus" : SurplusRefineProvider,
            "error" : ErrorRefineProvider,
            }

    eval = property(getEval, setEval, "Eval's Docstring")
    refine = property(getRefine, setRefine, "Refine's Docstring")