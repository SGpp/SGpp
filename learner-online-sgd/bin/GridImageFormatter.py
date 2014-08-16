##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
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

from bin.pysgpp import DataVector
from bin.learner.formatter import GridFormatter
import matplotlib.pyplot as plt 
from numpy import zeros, sqrt, ceil, floor


class GridImageFormatter(GridFormatter):
    '''
    classdocs
    '''


    ##Serializes grid to the stream
    #
    #@param memento: the Grid memento object
    #@param streamserializationStream: output stream where grid should be serialized to
    def serialize(self, memento, streamserializationStream):
        raise NotImplemented("The current method is not implemented in the class")
        
        
    ##Deserializes the Grid memento object from the stream
    #
    #@param serializationStream The stream to deserialize.
    #@return The Grid memento object.
    def deserialize(self, serializationStream):
        raise NotImplemented("The current method is not implemented in the class")
    
    
    ##Deserializes the Grid object from the file.
    # The file may or may be not gzip compressed.
    #@param filename The name of file with serialized Grid.
    #@return The Grid memento object.
    def deserializeFromFile(self, filename):
        raise NotImplemented("The current method is not implemented in the class")


    ##Serializes grid to the file
    #
    #@param memento: the Grid memento object
    #@param filename The name of file where the Grid object should be serialized to.
    #@todo (khakhutv) make the same variable names across similar functions
    def serializeToFile(self, memento, filename):
        fstream = self.gzOpen(filename, "w")

        try:
            grid = memento
            storage = grid.getStorage()
            coord_vector = DataVector(storage.dim())
            points = zeros([storage.size(), storage.dim()])
            for i in xrange(storage.size()):
                point = storage.get(i)
                point.getCoords(coord_vector)
                points[i] = [j for j in coord_vector.array()]

            rows = storage.dim()
            cols = storage.dim()

            i = 1

            for x1 in xrange(1,storage.dim()+1):
                 for x2 in xrange(1,storage.dim()+1):

                     plt.subplot(rows, cols, i)
                     plt.scatter(points[:,x1-1], points[:,x2-1])
                    
                     plt.axis([0,1,0,1])
                     plt.xticks([])
                     plt.yticks([])

                     if x1 == 1:
                         plt.title('x{}'.format(x2))

                     if x2 == 1:
                         plt.ylabel('x{}'.format(x1), rotation=90)

                     i +=1
            plt.tight_layout()
            plt.savefig(fstream)
            plt.close('all')
        finally:
            fstream.close()
