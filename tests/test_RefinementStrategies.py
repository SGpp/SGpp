###############################################################################
# $Copyright$ #
###############################################################################

import unittest



class TestRefinementANOVAStrategy(unittest.TestCase):
    def testRefine(self):
        from bin.pysgpp import HashRefinement, RefinementANOVAStrategy, Grid
        from bin.learner.LearnerBuilder import LearnerBuilder
        from bin.learner import Types
        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromARFFFile('refinement_strategy_sum')\
        .withGrid().withLevel(2)\
        .withBorder(Types.BorderTypes.TRAPEZOIDBOUNDARY)\
        .withSpecification().withIdentityOperator()\
        .withLambda(0.0001)\
        .withCGSolver().withImax(500)\
        .withStopPolicy().withAdaptiveItarationLimit(5)\
        .withProgressPresenter(InfoToScreen())\
        .andGetResult()
        
        learner.doLearningIteration(learner.dataContainer)
        
        
        
        
        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()