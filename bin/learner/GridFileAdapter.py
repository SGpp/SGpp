##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007-2009 Dirk Plueger (Dirk.Pflueger@in.tum.de)            #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)               #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

## @package GridFileAdapter
# @ingroup bin.learner
# @brief Grid File Adapter
# @version $CURR$

import re, gzip
from bin.pysgpp import *
from bin.learner import GridAdapter



class GridFileAdapter(GridAdapter):

    
    ##Load grid from the file
    #
    #@param filename: String with file name
    #@return: Grid
    def load(self, filename):
        fin = self.__gzOpen(filename, "r")
        text = fin.read()
        fin.close()
        return self.__restoreGrid(text)
    
    
    ##Unserializes grid from text
    # 
    #@param text: String with serialized grid information
    #@return: Grid object
    def __restoreGrid(self, text):
        #TODO: is there any control of correctness of text? (khakhutv)
        return Grid.unserialize(text)


    ##Save grid to the file
    #
    #@param grid: Grid object
    #@param destination: Filename the grid should be saved to
    def save(self, grid, destination):
        text = grid.serialize()
        fout = self.__gzOpen(destination, "w")
        fout.write(text)
        fout.close()
    
    
    ## Opens a file. If the file ends with ".gz", automatically gzip compression
    # is used for the file. Returns the filedescriptor
    # @param filename
    # @param mode, default: "r" for read only
    # @return file descriptor
    def __gzOpen(self, filename, mode="r"):
        # gzip-file
        if re.match(".*\.gz$", filename):
            # mode set for binary data?
            if not mode[-1] == "b":
                mode += "b"
            fd = gzip.open(filename, mode)
        # non gzip-file
        else:
            fd = open(filename, mode)
        return fd

