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


#from DataContainer import DataContainer

from time import time

class DataSpecification(object):
    """ generated source for DataSpecification

    """

    def isSaved(self):
        return self.__isSaved


    def setSaved(self):
        self.__isSaved = True

    filename = ''
    attributes = {}
    __isSaved = False
    
    def setFilename(self, filename):
        self.filename = filename
    
    def getFilename(self):
        #if filename is not set, it means the data subset was generated dynamically
        #so create new name, don't forget to save the data into the corresponding file!
        if self.filename == '':
            self.setFilename(self.generateFilename())
        return self.filename
    
    def generateFilename(self):
        return 'data.' + str(time()) + '.arff'
    
    def addAttribute(self, attributeName, attributeValue):
        self.attributes[attributeName] = attributeValue
        
    def getAttributes(self):
        return self.attributes
    
    def cleanAttributes(self):
        self.attributes = {}
        
    def createNumericAttributes(self, dim):
        self.cleanAttributes()
        for i in xrange(0,dim):
            self.addAttribute('x'+str(i), 'NUMERIC')
        self.addAttribute('class', 'NUMERIC')
    
    
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.
    # @todo (khakhutv) implement   
    def toString(self):
        serializationString = "'module' : '" + self.__module__ + "',\n"
        serializationString += '"filename":"' + self.getFilename() + '"\n'
        for key, value in self.getAttributes().items():
            serializationString = serializationString + ',"' + key + '" : "' + value + '"'
        return '{' + serializationString + '}'


