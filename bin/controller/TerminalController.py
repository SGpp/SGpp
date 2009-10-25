#############################################################################
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


import ConfigParser, os
from optparse import OptionParser

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.controller import InfoToScreen, InfoToFile, CheckpointController
from bin.learner.LearnerBuilder import LearnerBuilder
from bin.data.ARFFAdapter import ARFFAdapter


## The class provides the functionality for launching of learning jobs from
# terminal using console parameters or a configuration file.
#
# See file <i>tests/tbin/tcontroller/testsettings.job</i> for an example of configuration file.
# The command is then looks like:
# @code
# python TerminalController.py --jobfile path/to/file.job
# @endcode
# 
# In the same time there is an option to use command line for definition of
# all parameters of the job. Execute 
# @code
# python TerminalController.py --help
# @endcode
# for help to the console parameters.
#
class TerminalController:
 
    ## Initial processing of input parameters
    @classmethod
    def run(cls):
        # Initialize OptionParser, set Options
        parser = OptionParser()
        
        parser.add_option("--jobfile", action="store", type="string", dest="jobfile", help="Path to the file with job settings")
        
        parser.add_option("-l", "--level", action="store", type="int", dest="level", help="Gridlevel")
        parser.add_option("-a", "--adaptive", action="store", type="int", default=0, dest="adaptive", metavar="NUM", help="Using an adaptive Grid with NUM of refines")
        parser.add_option("--adapt_points", action="store", type="int", default=1, dest="adapt_points", metavar="NUM", help="Number of points in one refinement iteration")
        parser.add_option("--adapt_rate", action="store", type="float", dest="adapt_rate", metavar="NUM", help="Percentage of points from all refinable points in one refinement iteration")
        parser.add_option("--adapt_start", action="store", type="int", default=0, dest="adapt_start", metavar="NUM", help="The index of adapt step to begin with")
        parser.add_option("--adapt_threshold", action="store", type="float", default=0.0, dest="adapt_threshold", metavar="NUM", help="The threshold, an error or alpha has to be greater than in order to be reined.")
        parser.add_option("-m", "--mode", action="store", type="string", default="apply", dest="mode", help="Specifies the action to do. Get help for the mode please type --mode help.")
        parser.add_option("-C", "--zeh", action="store", type="string", default="laplace", dest="zeh", help="Specifies the action to do.")
        parser.add_option("-f", "--foldlevel", action="store", type="int",default=10, metavar="LEVEL", dest="f_level", help="If a fold mode is selected, this specifies the number of sets generated")
        parser.add_option("--onlyfoldnum", action="store", type="int", default=-1, metavar="I", dest="onlyfoldnum", help="Run only fold I in n-fold cross-validation. Default: run all")
        parser.add_option("-L", "--lambda", action="store", type="float",default=0.000001, metavar="LAMBDA", dest="regparam", help="Lambda")
        parser.add_option("-i", "--imax", action="store", type="int",default=500, metavar="MAX", dest="imax", help="Max number of iterations")
        parser.add_option("-r", "--accuracy", action="store", type="float",default=0.0001, metavar="ACCURACY", dest="r", help="Specifies the accuracy of the CG-Iteration")
        parser.add_option("--max_accuracy", action="store", type="float", default=None, metavar="ACCURACY", dest="max_r", help="If the norm of the residuum falls below ACCURACY, stop the CG iterations")
        parser.add_option("-d", "--data", action="append", type="string", dest="data", help="Filename for the Datafile.")
        parser.add_option("-t", "--test", action="store", type="string", dest="test", help="File containing the testdata")
        parser.add_option("-A", "--alpha", action="store", type="string", dest="alpha", help="Filename for a file containing an alpha-Vector")
        parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="Filename where the calculated alphas are stored")
        parser.add_option("-g", "--gnuplot", action="store", type="string", dest="gnuplot", help="In 2D case, the generated can be stored in a gnuplot readable format.")
        parser.add_option("-R", "--resolution", action="store", type="int",default=50, metavar="RESOLUTION", dest="res", help="Specifies the resolution of the gnuplotfile")
        parser.add_option("-s", "--stats", action="store", type="string", dest="stats", help="In this file the statistics from the test are stored")
        parser.add_option("-p", "--polynom", action="store", type="int", default=0, dest="polynom", help="Sets the maximum degree for high order basis functions. Set to 2 or larger to activate. Works only with 'identity' and 'fold'-modes.")
        parser.add_option("-b", "--border", action="store_true", default=False, dest="border", help="Enables special border base functions")
        parser.add_option("--trapezoid-boundary", action="store_true", default=False, dest="trapezoidboundary", help="Enables boundary functions that have a point on the boundary for every inner point (Trapezoid)")
        parser.add_option("--complete-boundary", action="store_true", default=False, dest="completeboundary", help="Enables boundary functions that have more points on the boundary than inner points")
        parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="Provides extra output")
        parser.add_option("--normfile", action="store", type="string", dest="normfile", metavar="FILE", help="For all modes that read data via stdin. Normalizes data according to boundaries in FILE")
        parser.add_option("--reuse", action="store_true", default=False, dest="reuse", help="Reuse alpha-values for CG")
        parser.add_option("--seed", action="store", type="float", dest="seed", help="Random seed used for initializing")
        parser.add_option("--regression", action="store_true", default=False, dest="regression", help="Use regression approach.")
        parser.add_option("--checkpoint", action="store", type="string", dest="checkpoint", help="Filename for checkpointing. For fold? and test. No file extension.")
        parser.add_option("--grid", action="store", type="string", dest="grid", help="Filename for Grid-resume. For fold? and test. Full filename.")
        parser.add_option("--epochs_limit", action="store", type="int", default="0", dest="epochs_limit", help="Number of refinement iterations (epochs), MSE of test data have to increase, before refinement will stop.")
        parser.add_option("--mse_limit", action="store", type="float", default="0.0", dest="mse_limit", help="If MSE of test data fall below this limit, refinement will stop.")
        parser.add_option("--grid_limit", action="store", type="int", default="0", dest="grid_limit", help="If the number of points on grid exceed grid_limit, refinement will stop.")
    
        # parse options
        (options,args)=parser.parse_args()
        
        #if job defined in file:
        if(options.jobfile != None):
            TerminalController.constructObjectsFromFile(options.jobfile)
        #if job defined with arguments:
        else:
            TerminalController.constructObjectsFromOptions(options)
        
    
    ## Construct all objects from the information, defined in the job file
    #
    # @param filename: string path to the file with job settings
    @classmethod
    def constructObjectsFromFile(cls, filename):
        configuration = ConfigParser.ConfigParser()
        configuration.readfp(open(filename, 'r'))
        
        #Create builder for specified learner
        builder = LearnerBuilder()
        learner_type = configuration.get('learner', 'type')
        if(learner_type == 'classification'):
            builder.buildClassifier()
        elif(learner_type == 'regression'):
            builder.buildRegressor()
        else:
            raise Exception('Wrong learner type in job configuration file')
        
        
        #dataset options
        options = TerminalController.itemsToDict(configuration.items('data'))
        
        if options['file_type'] == 'arff':
            if options.has_key('train_file'): builder.withTrainingDataFromARFFFile(options['train_file'])
            else: raise Exception('Path to file with training data set is not defined in configurationfile')
            
            if options.has_key('test_file'): builder.withTestingDataFromARFFFile(options['test_file'])
       
        else: raise Exception('Unsupported data type in job configuration file')
        
        #grid options
        builder = builder.withGrid()
        options = TerminalController.itemsToDict(configuration.items('grid'))
        
        if options.has_key('grid_file'):
            builder.fromFile(options['grid_file'])
        else:
            try:
                if options.has_key('level'): builder.withLevel(int(options['level']))
                if options.has_key('polynomial'): builder.withPolynomialBase(int(options['polynomial']))
                if options.has_key('border'): builder.withBorder(options['border'])
            except: raise Exception('Grid configuration in job file is incorrect')
            
        #training specification
        builder = builder.withSpecification()
        options = TerminalController.itemsToDict(configuration.items('refinement'))
        
        if options.has_key('points'): 
            if '.' in options['points']:
                builder.withAdaptRate( float(options['points']) )
            else:
                builder.withAdaptPoints( int(options['points']) )
                
        options = TerminalController.itemsToDict(configuration.items('learner'))
        
        if options.has_key('regularization_parameter'): builder.withLambda( float(options['regularization_parameter']) )
        if options.has_key('regularization_operator'):
            if options['regularization_operator'] == 'laplace': builder.withLaplaceOperator()
            elif options['regularization_operator'] == 'idenitty': builder.withIdentityOperator()
            else: raise Exception('Incorrect regulariation operator type')
        
        if options.has_key('threshold'): builder.withAdaptThreshold( float(options['threshold']) )
        
        
        #stop policy
        builder = builder.withStopPolicy()
        options = TerminalController.itemsToDict(configuration.items('refinement'))
        
        if options.has_key('iterations'): builder.withAdaptiveItarationLimit( int(options['iterations']) )
        if options.has_key('gridsize'): builder.withGridSizeLimit( int(options['gridsize']) )
        if options.has_key('mse'): builder.withMSELimit( float(options['mse']) )
        if options.has_key('epochs'): builder.withEpochsLimit( int(options['epochs']) )
        
        # linear solver
        builder = builder.withCGSolver()
        options = TerminalController.itemsToDict(configuration.items('solver'))
        
        if options.has_key('accuracy'): builder.withAccuracy( float(options['accuracy']) )
        if options.has_key('imax'): builder.withImax( int(options['imax']) )
        
        #presentor
        options = TerminalController.itemsToDict(configuration.items('output'))
        if options.has_key('type'):
            types = options['type'].split(',')
            for type in types:
                if type.strip() == 'InfoToScreen': builder.withProgressPresentor(InfoToScreen())
                if type.strip() == 'InfoToFile':
                    if options.has_key('filename'): builder.withProgressPresentor(InfoToFile(options['filename']))
                    else: raise Exception('Define filename in order to use InfoToFile output')
                    
        #checkpoint
        options = TerminalController.itemsToDict(configuration.items('checkpoints'))
        
        if options.has_key('name'):
            checkpointController = CheckpointController(options['name'], options['path'], int(options['interval']) )
                                                                                                      
            builder.withCheckpointController(checkpointController)
            
        # @todo (khakhutv) develop an intelligent way to restore from checkpoints
            
        #Get Results and perform wanted action
        learner = builder.andGetResult()
        
        options = TerminalController.itemsToDict(configuration.items('learner'))
        
        if options['action'] == 'learn':
            if not options['with_testing'] or options['with_testing'].lower() == 'no':
                learner.learnData()
            elif options['with_testing'].lower() == 'yes': learner.learnDataWithTest()
            else: raise Exception('with_testion can only be set to "yes" or "no"')
        elif options['action'] == 'apply':
            points_file = configuration.get('data', 'points_file')
            if points_file != None:
                learner.applyData( ARFFAdapter(points_file).loadData().getPoints() )
            else: raise Exception('To evaluate value of points define file path "points_file" in the section "data"')
        else: raise Exception('Incorrect action in job configuration file')
    
    
    
    
    ## Construct all objects form the information, defined with arguments
    #
    # @param options: OptionParser result object with options
    @classmethod
    def constructObjectsFromOptions(cls, options):
        #@todo (khakhutv) implement terminal initialization from command line parameters
        return 
    
    
    
    
    
    ## Convert list of items [(key1,value1), ...] to dictionary
    #
    # @param items: list of items
    # @return: dictionary
    @classmethod
    def itemsToDict(cls, items):
        result = {}
        for item in items:
            result[item[0]] = item[1]
            
        return result
    
    

if __name__=='__main__':
    TerminalController.run()
