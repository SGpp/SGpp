# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#############################################################################
                                    #
#############################################################################
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer


import ConfigParser, os
from optparse import OptionParser

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.controller import InfoToScreen, InfoToScreenRegressor, InfoToFile, CheckpointController
from pysgpp.extensions.datadriven.learner.LearnerBuilder import LearnerBuilder
from pysgpp.extensions.datadriven.data.ARFFAdapter import ARFFAdapter
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
import types


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
        parser.add_option("--generate", action="store_true", default=False, dest="generate", help="Only generate the code for LearnerBuilder, but don't evaluate it")

        # parse options
        (options,args)=parser.parse_args()
        
        #if job defined in file:
        if(options.jobfile != None):
            TerminalController.constructObjectsFromFile(options.jobfile)
        #if job defined with arguments:
        elif options.generate:
            print TerminalController.generateBuilderCodeFromOptions(options)
        else:
            TerminalController.constructObjectsFromOptions(options)
        
    
    ## Construct all objects from the information, defined in the job file
    #
    # @param cls python keyword for class method (no specification needed)
    # @param filename: string path to the file with job settings
    @classmethod
    def constructObjectsFromFile(cls, filename):
        configuration = ConfigParser.ConfigParser()
        configuration.readfp(open(filename, 'r'))
        
        #Create builder for specified learner
        builder = LearnerBuilder()
        learner_type = configuration.get('learner', 'type')
        if learner_type == 'classification':
            builder.buildClassifier()
        elif learner_type == 'regression':
            builder.buildRegressor()
        else:
            raise Exception('Wrong learner type in job configuration file')
        
        
        #dataset options
        options = TerminalController.itemsToDict(configuration.items('data'))
        
        if options['file_type'] == 'arff':
            if options.has_key('train_file'): 
                if type(options['train_file']) != types.ListType:
                    builder.withTrainingDataFromARFFFile(options['train_file'])
                else:
                    fileCounter = 0
                    for train_file in options['train_file']:
                        builder.withTrainingDataFromARFFFile(train_file, DataContainer.TRAIN_CATEGORY + str(fileCounter))
                        fileCounter += 1
                                                            
                        
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
        if options.has_key('max_threshold'): builder.withThreshold(float(options['max_threshold']))
        
        #presentor
        options = TerminalController.itemsToDict(configuration.items('output'))
        if options.has_key('type'):
            types = options['type'].split(',')
            for type in types:
                if type.strip() == 'InfoToScreen': builder.withProgressPresenter(InfoToScreen())
                if type.strip() == 'InfoToScreenRegressor': builder.withProgressPresenter(InfoToScreenRegressor())
                elif type.strip() == 'InfoToFile':
                    if options.has_key('filename'): builder.withProgressPresenter(InfoToFile(options['filename']))
                    else: raise Exception('Define filename in order to use InfoToFile output')
                    
        #checkpoint
        options = TerminalController.itemsToDict(configuration.items('checkpoints'))
        
        if options.has_key('name'):
            path = options['path'] if options.has_key('path') else None
            interval = options['interval'] if options.has_key('inte    rval') else None
            checkpointController = CheckpointController(options['name'], path, interval)
            if options.has_key('restore_iteration'):
                learner = checkpointController.loadAll(options['restore_iteration'])
                                                                                                      
            builder.withCheckpointController(checkpointController)
            
            
        #Get Results and perform wanted action
        # if learner was not created by checkpoint controller, create it with learner builder
        try:
            if learner not in dir():
                learner = builder.andGetResult()
        except: learner = builder.andGetResult()
        
        # Folding
        options = TerminalController.itemsToDict(configuration.items('folding'))
        if options.type in ['fold', 'folds', 'foldstratified', 'foldr', 'foldf'] :
            if options.type == 'fold': builder.withRandomFoldingPolicy() 
            elif options.type == 'folds': builder.withSequentialFoldingPolicy()
            elif options.type in ['foldstratified', 'foldr']: builder.withStratifiedFoldingPolicy()
            elif options.type == 'foldf': builder.withFilesFoldingPolicy()
            if options.seed: builder.withSeed(options.seed)
            if options.level: builder.withLevel(options.level)
        else: raise Exception('Unknown folding type')
            
        
        options = TerminalController.itemsToDict(configuration.items('learner'))
        if options['action'] == 'learn':
            if (not options.has_key('with_testing')) or options['with_testing'].lower() == 'no':
                learner.learnData()
            elif options['with_testing'].lower() == 'yes': learner.learnDataWithTest()
            else: raise Exception('with_testion can only be set to "yes" or "no"')
        elif options['action'] == 'apply':
            points_file = configuration.get('data', 'points_file')
            if points_file != None:
                learner.applyData( ARFFAdapter(points_file).loadData().getPoints() )
            else: raise Exception('To evaluate value of points define file path "points_file" in the section "data"')
        elif options['action'] == 'fold':
            builder.getCheckpointController().generateFoldValidationJob('PUT_YOUR_EMAIL_HERE')
        else: raise Exception('Incorrect action in job configuration file')
    
    
    
    
    ## Construct all objects form the information, defined with arguments
    #
    # @param cls python keyword for class method (no specification needed)
    # @param options: OptionParser result object with options
    @classmethod
    def constructObjectsFromOptions(cls, options):       
        #Create builder for specified learner
        builder = LearnerBuilder()
        if  options.regression:
            builder.buildRegressor()

        else:
            builder.buildClassifier()
        
        # load alpha file
        if options.alpha:
            builder.withInitialAlphaFromARFFFile(options.alpha)
        
        
        #dataset options
        if len(options.data) == 1:
            builder.withTrainingDataFromARFFFile(options.data[0])
        elif len(options.data) > 1:
            fileCounter = 0
            for filename in options.data:
                builder.withTrainingDataFromARFFFile(filename, DataContainer.TRAIN_CATEGORY + str(fileCounter))
                fileCounter += 1
        else: 
            raise Exception('Define the path to the training data set')
            
        if options.test: builder.withTestingDataFromARFFFile(options.test)
       
        
        #grid options
        builder = builder.withGrid()
        if options.grid:
            builder.fromFile(options.grid)
        else:
            try:
                if options.level: builder.withLevel(options.level)
                if options.polynom: builder.withPolynomialBase(options.polynom)
                if options.trapezoidboundary: builder.withBorder(BorderTypes.TRAPEZOIDBOUNDARY)
                elif options.completeboundary: builder.withBorder(BorderTypes.COMPLETEBOUNDARY)
                # @fixme (khakhutv)the name "NONE" for the border type should be changed to something more meaningful
                elif options.border: builder.withBorder(BorderTypes.NONE)
            except: raise Exception('Grid configuration arguments incorrect')
        
        if options.adapt_start: builder.withStartingIterationNumber(options.adapt_start)
            
        #training specification
        builder = builder.withSpecification()
        if options.adapt_rate: builder.withAdaptRate(options.adapt_rate)
        elif options.adapt_points: builder.withAdaptPoints(options.adapt_points)
                
        if options.regparam: builder.withLambda(options.regparam)
        if options.zeh:
            if options.zeh == 'laplace': builder.withLaplaceOperator()
            elif options.zeh == 'identity': builder.withIdentityOperator()
            else: raise Exception('Incorrect regulariation operator type')
        
        if options.adapt_threshold: builder.withAdaptThreshold(options.adapt_threshold)
        
        
        #stop policy
        builder = builder.withStopPolicy()
        if options.adaptive: builder.withAdaptiveItarationLimit(options.adaptive)
        if options.grid_limit: builder.withGridSizeLimit(options.grid_limit)
        if options.mse_limit: builder.withMSELimit(options.mse_limi)
        if options.epochs_limit: builder.withEpochsLimit(options.epochs_limit)
        
        # linear solver
        builder = builder.withCGSolver()        
        if options.r: builder.withAccuracy(options.r)
        if options.max_r: builder.withThreshold(options.max_r)
        if options.imax: builder.withImax(options.imax)
        
        #presentor
        if options.verbose: # print to the screen
            if options.regression:
                builder.withProgressPresenter(InfoToScreenRegressor())
            else:
                builder.withProgressPresenter(InfoToScreen())
        if options.stats:
            builder.withProgressPresenter(InfoToFile(options.stats))
            
                    
        #checkpoint
        if options.checkpoint:
            title = os.path.basename(options.checkpoint)
            path = os.path.dirname(options.checkpoint)
            checkpointController = CheckpointController(title, path)
                                                                                                      
            builder.withCheckpointController(checkpointController)
            
        # Folding
        if options.mode in ['fold', 'folds', 'foldstratified', 'foldf', 'foldr'] :
            if options.mode == 'fold': builder.withRandomFoldingPolicy()
            elif options.mode == 'folds': builder.withSequentialFoldingPolicy()
            elif options.mode in ['foldstratified', 'foldr']: builder.withStratifiedFoldingPolicy()
            elif options.mode == 'foldf': builder.withFilesFoldingPolicy()
            if options.seed: builder.withSeed(options.seed)
            if options.level: builder.withLevel(options.level)

            
            
        #Get Results and perform wanted action
        learner = builder.andGetResult()
        options.mode = options.mode.lower()
        if options.mode == 'normal':
            learner.learnData()
        elif options.mode == 'test':
            learner.learnDataWithTest()
        elif options.mode == 'apply':
            learner.applyData( ARFFAdapter(options.data).loadData().getPoints() )
        elif options.mode in ['fold', 'folds', 'foldstratified', 'foldf', 'foldr']:
            builder.getCheckpointController().generateFoldValidationJob('PUT_YOUR_EMAIL_HERE')
        elif options.mode in ['eval', 'evalstdin']:
            raise Exception('This action is not implemented yet')
        else: raise Exception('Incorrect action configuration')
    
    
    
    
    
    ## Convert list of items [(key1,value1), ...] to dictionary
    #
    # @param cls python keyword for class method (no specification needed)
    # @param items: list of items
    # @return: dictionary
    @classmethod
    def itemsToDict(cls, items):
        result = {}
        for item in items:
            if result.has_key(item[0]) == False:
                result[item[0]] = item[1]
            else:
                # the old value is already a list
                if type( result[item[0]] ) == types.ListType:
                    result[ item[0] ].append( item[1] )
                # create new list with old and new values
                else:
                    result[item[0]] = [ result[item[0]], item[1] ]
            
            
        return result
    
    

    ## Generate LearnerBuilder code according to parameters
    # You can execute the code with exec command and use builder object in your 
    # python code
    #
    # @param cls python keyword for class method (no specification needed)
    # @param options: OptionParser result object with options
    @classmethod
    def generateBuilderCodeFromOptions(cls, options):    
        #Create builder for specified learner
        code = "builder = LearnerBuilder()\n"
        if  options.regression:
            code += "builder.buildRegressor()\n"

        else:
            code += "builder.buildClassifier()\n"
            
        # load alpha file
        if options.alpha:
            code += "builder.withInitialAlphaFromARFFFile('%s')\n" % options.alpha
        
        
        #dataset options
        if len(options.data) == 1:
            code += "builder.withTrainingDataFromARFFFile('%s')\n" % options.data[0]
        elif len(options.data) > 1:
            fileCounter = 0
            for filename in options.data:
                code += "builder.withTrainingDataFromARFFFile('%s', '%s')\n" % (filename, DataContainer.TRAIN_CATEGORY + str(fileCounter))
                fileCounter += 1
        else: 
            raise Exception('Define the path to the training data set')
            
        if options.test: code += "builder.withTestingDataFromARFFFile('%s')\n" % options.test
       
        
        #grid options
        code += "builder = builder.withGrid()\n"
        if options.grid:
            code += "builder.fromFile('%s')\n" % options.grid
        else:
            try:
                if options.level: code += "builder.withLevel(%d)\n" % options.level
                if options.polynom: code += "builder.withPolynomialBase(%d)\n" % options.polynom
                if options.trapezoidboundary: code += "builder.withBorder(BorderTypes.TRAPEZOIDBOUNDARY)\n"
                elif options.completeboundary: code += "builder.withBorder(BorderTypes.COMPLETEBOUNDARY)\n"
                # @fixme (khakhutv)the name "NONE" for the border type should be changed to something more meaningful
                elif options.border: code += "builder.withBorder(BorderTypes.NONE)\n"
            except: raise Exception('Grid configuration arguments incorrect')
            
            
        if options.adapt_start: code += "builder.withStartingIterationNumber(%d)\n" % options.adapt_start 
            
        #training specification
        code += "builder = builder.withSpecification()\n"
        if options.adapt_rate: code += "builder.withAdaptRate(%d)\n" % options.adapt_rate
        elif options.adapt_points: code += "builder.withAdaptPoints(%d)\n" % options.adapt_points
                
        if options.regparam: code += "builder.withLambda(%1.20f)\n" % options.regparam
        if options.zeh:
            if options.zeh == 'laplace': code += "builder.withLaplaceOperator()\n"
            elif options.zeh == 'identity': code += "builder.withIdentityOperator()\n"
            else: raise Exception('Incorrect regulariation operator type')
        
        if options.adapt_threshold: code += "builder.withAdaptThreshold(%f)\n" % options.adapt_threshold
        
        
        #stop policy
        code += "builder = builder.withStopPolicy()\n"
        if options.adaptive: code += "builder.withAdaptiveItarationLimit(%d)\n" % options.adaptive
        if options.grid_limit: code += "builder.withGridSizeLimit(%d)\n" % options.grid_limit
        if options.mse_limit: code += "builder.withMSELimit(%f)\n" % options.mse_limit
        if options.epochs_limit: code += "builder.withEpochsLimit(%d)\n" % options.epochs_limit
        
        # linear solver
        code += "builder = builder.withCGSolver()  \n"      
        if options.r: code += "builder.withAccuracy(%f)\n" % options.r
        if options.max_r: code += "builder.withThreshold(%f)\n" % options.max_r
        if options.imax: code += "builder.withImax(%d)\n" % options.imax
        
        #presentor
        if options.verbose: # print to the screen
            if options.regression:
                code += "builder.withProgressPresenter(InfoToScreenRegressor())\n"
            else:
                code += "builder.withProgressPresenter(InfoToScreen())\n"
        if options.stats:
            code += "builder.withProgressPresenter(InfoToFile('%s'))\n" % options.stats
            
                    
        #checkpoint
        if options.checkpoint:
            title = os.path.basename(options.checkpoint)
            path = os.path.dirname(options.checkpoint)
            code += "checkpointController = CheckpointController('%s', '%s')\n" % (title, path)
            code += "builder.withCheckpointController(checkpointController)\n"
            
        # Folding
        if options.mode in ['fold', 'folds', 'foldstratified', 'foldf', 'foldr'] :
            if options.mode == 'fold': code += "builder.withRandomFoldingPolicy()\n"
            elif options.mode == 'folds': code += "builder.withSequentialFoldingPolicy()\n"
            elif options.mode in ['foldstratified', 'foldr']: code += "builder.withStratifiedFoldingPolicy()\n"
            elif options.mode == 'foldf': code += "builder.withFilesFoldingPolicy()\n"
            if options.seed: code += "builder.withSeed(%d)\n" % options.seed
            if options.level: code += "builder.withLevel(%d)\n" % options.level

            
        return code
        
        
        
        
if __name__=='__main__':
    TerminalController.run()