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
from bin.tools import *

class EvalProvider:
    
    def __init__(self, status):
        self.status = status
        
        self.training_overall = []
        self.testing_overall = []
        self.number_points = []
        
    def reset(self):
        """Called before each learning step"""
        self.training_results = []
        self.testing_results = []
        
    def updateResults(self, alpha, trainingData, testingData):
        """Calculate new results for training and test data for statistics"""
        self.training(alpha, trainingData)
        self.testing(alpha, testingData)
        
class ClassesEvalProvider(EvalProvider):
    
    def testVector(self, gridDecorator, alpha, data, classes):
        testOP = gridDecorator.grid.createOperationTest()
        return testOP.test(alpha, data, classes) / float(data.getSize())
        
    def training(self, alpha, data):
        """Called to evaluate training data once per fold. Used later for statistics."""
        if data:
            self.c.append(self.testVector(self.status, alpha, data[0], data[1]))

    def testing(self, alpha, data):
        """Called to evaluate testing data once per fold. Used later for statistics."""
        if data:
            self.testing_results.append(self.testVector(self.status, alpha, data[0], data[1]))
                                     
    def updateStatistics(self,  adaptive, l, level, stats, verbose):
        """Called after each learning step"""
        i = float(len(self.training_results))
        self.training_overall.append(sum(self.training_results)/i)
        self.testing_overall.append(sum(self.testing_results)/i)
        self.number_points.append(self.status.getSize())
        
        if verbose:       
            print "Traing avg MSE:"
            print self.training_overall
            print "Testing avg. MSE:"
            print self.testing_overall

        if stats != None:
            txt = "%f, %-10g, %f" % (level, l, adaptive)
            for i in xrange(len(self.training_overall)):
                txt = txt + ", %f, %.10f, %.10f" % (self.number_points[i], self.training_overall[i], self.testing_overall[i])
            if verbose:
                print "Statistics: level, lambda, adaptive, number of point, training overall, testing overall"
                print txt
            writeLockFile(stats, txt+"\n")


class RegressionEvalProvider(EvalProvider):

    def training(self, alpha, data):
        """Called to evaluate training data once per fold. Used later for statistics."""
        if data:
            self.training_results.append(self.regressionTest(alpha, data[0], data[1]))

    def testing(self, alpha, data):
        """Called to evaluate testing data once per fold. Used later for statistics."""
        if data:
            self.testing_results.append(self.regressionTest(alpha, data[0], data[1]))
             
    def regressionTest(self, alpha, data, classes):
        temp = DataVector(classes.getSize())
        B = self.status.grid.createOperationB()
        B.multTranspose(alpha, data, temp)
        
        temp.sub(classes)
        temp.sqr()
        return temp.sum()/temp.getSize()
        
                                     
    def updateStatistics(self,  adaptive, l, level, stats, verbose):
        """Called after each learning step"""
        i = float(len(self.training_results))
        self.training_overall.append(sum(self.training_results)/i)
        self.testing_overall.append(sum(self.testing_results)/i)
        self.number_points.append(self.status.getSize())
        
        if verbose:       
            print "Traing avg MSE:"
            print self.training_overall
            print "Testing avg. MSE:"
            print self.testing_overall
        
#        if options.verbose:
#            print "MSE: %2.10f" %self.testing_results[-1]

        if stats != None:
            txt = "%f, %-10g, %f" % (level, l, adaptive)
            for i in xrange(len(self.training_overall)):
                txt = txt + ", %f, %.10f, %.10f" % (self.number_points[i], self.training_overall[i], self.testing_overall[i])
            if verbose:
                print "Statistics: level, lambda, adaptive, number of point, training overall, testing overall"
                print txt
            writeLockFile(stats, txt+"\n")