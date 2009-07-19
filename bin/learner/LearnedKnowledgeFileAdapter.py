##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2008 Dirk Plueger (pflueged@in.tum.de)                      #
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

## @package LearnedKnowledgeFileAdapter
# @ingroup bin.learner
# @brief Storing and restoring learned knowledge in files
# @version $CURR$

import sys, re, gzip

from LearnedKnowledge import LearnedKnowledge
from KnowledgeAdapter import KnowledgeAdapter
from bin.pysgpp import DataVector

class LearnedKnowledgeFileAdapter(KnowledgeAdapter):


    ## Load knowledge data from the source
    #
    # @param filename: String with filename it should be loaded from
    def load(self, filename):
        if type(filename) != type(''): raise AttributeError, "Filename as destination expected"
        alphas = self.__readAlphaARFF(filename)
        knowledge = LearnedKnowledge(self)
        knowledge.update(alphas)
        return knowledge
        
        
    ## Save knowledge data to the file
    #
    # @param knowledge: LearnedKnowledge object to save
    # @param dest: String file name
    # @return: LearnedKnowledge object
    def save(self, knowledge, dest):
        if type(dest) != type(''): raise AttributeError, "Filename as destination expected"
        alphas = knowledge.getAlphas();
        self.__writeAlphaARFF(dest, alphas)
        return knowledge
    
    
    ## Reads alpha vector from the ARFF file
    #
    # @param filename: String file name
    # @param alpha: DataVector alpha file
    # @return: alpha DataVector
    def __readAlphaARFF(self, filename):
        data = self.__readDataARFF(filename)
        alpha = DataVector(len(data["data"][0]), 1)
        for i in xrange(len(data["data"][0])):
            alpha[i] = data["data"][0][i]
        return alpha
    
    
    ## Writes alpha vector to the ARFF file
    #
    # @param filename: String file name
    # @param alpha: DataVector of alpha
    def __writeAlphaARFF(self, filename, alpha):
        fout = self.__gzOpen(filename, "w")
        fout.write("@RELATION \"%s ALPHAFILE\"\n\n" % filename)
        fout.write("@ATTRIBUTE alpha NUMERIC\n")
        
        fout.write("\n@DATA\n")
        
        for i in xrange(len(alpha)):
            fout.write("%1.20f\n" % alpha[i])
        
        fout.close()
        
        
    ## Reads in an ARFF file:
    # The data is stored in lists. There is a value list for every dimension of the data set. e.g. 
    # [[2, 3],[1, 1]] are the data points P_1(2,1) and P_2(3,1)
    #
    # @param filename the file's filename that should be read
    # @return returns a set of a array with the data (named data), a array with the classes (named classes) and the filename named as filename
    def __readDataARFF(self, filename):
        fin = self.__gzOpen(filename, "r")
        data = []
        classes = []
        hasclass = False
    
        # get the different section of ARFF-File
        for line in fin:
            sline = line.strip().lower()
            if sline.startswith("%") or len(sline) == 0:
                continue
    
            if sline.startswith("@data"):
                break
            
            if sline.startswith("@attribute"):
                value = sline.split()
                if value[1].startswith("class"):
                    hasclass = True
                else:
                    data.append([])
        
        #read in the data stored in the ARFF file
        for line in fin:
            sline = line.strip()
            if sline.startswith("%") or len(sline) == 0:
                continue
    
            values = sline.split(",")
            if hasclass:
                classes.append(float(values[-1]))
                values = values[:-1]
            for i in xrange(len(values)):
                data[i].append(float(values[i]))
                
        # cleaning up and return
        fin.close()
        return {"data":data, "classes":classes, "filename":filename}
        
        
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

