
'''
Created on Feb 6, 2015

@author: franzefn
'''

from interpolationAlgorithm import InterpolationAlgorithm
from bin.uq.dists.SGDEdist import SGDEdist
import ConfigParser as cp


class EstimateDensityAlgorithm(InterpolationAlgorithm):

    def __init__(self, config):
        # load path and file names of the grid and if the coefficients from
        # the config file
        s = cp.ConfigParser()
        s.optionxform = str  # force the parser to read keywords case sensitive
        s.read(config)
        if 'gridFile' not in s.options('dmest'):
            raise AttributeError("option 'gridFile' is required in section 'dmest'")
        if 'writeGridFile' not in s.options('dmest'):
            raise AttributeError("option 'writeGridFile' is required in section 'dmest'")
        if 'writeAlphaFile' not in s.options('dmest'):
            raise AttributeError("option 'writeAlphaFile' is required in section 'dmest'")

        # load the file where the grid is stored
        self.gridfile = s.get('dmest', 'writeGridFile')

        if self.gridfile != s.get('dmest', 'gridFile'):
            raise AttributeError("option 'writeGridFile' has to be the same as 'gridFile'")

        # remove refinement if there is some
        if "refNr" in s.options("dmest"):
            s.set('dmest', 'refNr', 0)

        # write new config to file
        tmpConfig = "%s.tmp" % config
        fd = open(tmpConfig, "w")
        s.write(fd)
        fd.close()

        self.config = tmpConfig

    def computeHierarchicalCoefficients(self, grid, alpha, newGridPoints):
        # serialize grid
        fd = open(self.gridfile, "w")
        fd.write(grid.serialize())
        fd.close()

        # estimate the density on the given grid
        dist = SGDEdist.byConfig(self.config)

        # return new new coefficient vector
        return dist.alpha
