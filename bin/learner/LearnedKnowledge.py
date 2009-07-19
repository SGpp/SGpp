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

## @package LearnedKnowledge
# @ingroup bin.learner
# @brief Represents learned knowledge 
# @version $CURR$


class LearnedKnowledge(object):

    
    __alphas = None         #DataVector with learned alpha
    __dataSource = None     #KnowledgeAdapter, which handles loading and restoring
    
    
    ##Constructor
    #
    #@param adapter: KnowledgeAdapter, default value is None
    def __init__(self, adapter=None):
        self.__dataSource = adapter


    ##Save object to the destination
    #
    #@param destination: Destination, i.e. filename
    def save(self, destintation):
        return self.__dataSource.save(self, destintation)
    
    
    ##Load LearnedKnowledge from the source
    #
    #@param source: KnowledgeAdapter object
    #@return: Restored LearnedKnowledge object
    @classmethod
    def load(cls, source):
        knowledge = source.load()
        return knowledge


    ##Returns current alpha vector
    #
    #@return: DataVector of current alpha
    def getAlphas(self):
        return self.__alphas


    ##Alters the current alpha vector
    #
    #@param alphas: new alpha vector
    #@return: LearnedKnowledge itself
    def update(self, alpha):
        self.__alphas = alpha
        return self


