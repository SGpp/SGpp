/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/application/Learner.hpp"

namespace sg
{

namespace datadriven
{

Learner::Learner(bool verbose) : alpha_(NULL), grid_(NULL), verbose_(verbose)
{
}

Learner::~Learner()
{
	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;
}

void Learner::createInitialGrid(sg::base::RegularGridConfiguration& GridConfig)
{
	if (GridConfig.type_ == sg::base::LinearTrapezoidBoundary)
	{
		grid_ = new sg::base::LinearTrapezoidBoundaryGrid(GridConfig.dim_);
	}
	else if (GridConfig.type_ == sg::base::ModLinear)
	{
		grid_ = new sg::base::ModLinearGrid(GridConfig.dim_);
	}
	else if (GridConfig.type_ == sg::base::Linear)
	{
		grid_ = new sg::base::LinearGrid(GridConfig.dim_);
	}
	else
	{
		std::cout << std::endl << "An unsupported grid type was chosen! Exiting...." << std::endl << std::endl;
		grid_ = NULL;
		return;
	}

	// Generate regular Grid with LEVELS Levels
	sg::base::GridGenerator* myGenerator = grid_->createGridGenerator();
	myGenerator->regular(GridConfig.level_);
	delete myGenerator;
}

void Learner::trainGrid(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes, sg::base::DataMatrix& testDataset,
		sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal,
		sg::base::AdpativityConfiguration& AdaptConfig, sg::base::OperationMatrix& SLESystem)
{
    // Generate CG to solve System
	sg::solver::SLESolver* myCG;

	if (SolverConfigRefine.type_ != sg::solver::CG)
	{
		myCG = new sg::solver::ConjugateGradients(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
	}
	else if (SolverConfigRefine.type_ != sg::solver::BiCGSTAB)
	{
		myCG = new sg::solver::BiCGStab(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
	}
	else
	{
		// @TODO: Error
	}

    std::cout << "Starting Learning...." << std::endl;

    // execute adaptsteps
    sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();

    for (size_t i = 0; i < AdaptConfig.noRefinements_+1; i++)
    {
    	std::cout << std::endl << "Doing refinement: " << i << std::endl;

    	myStopwatch->start();

    	// Do Refinements
    	if (i > 0)
    	{
    		sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(alpha_, AdaptConfig.noPoints_, AdaptConfig.threshold_);
    		grid_->createGridGenerator()->refine(myRefineFunc);
    		delete myRefineFunc;

    		SLESystem.rebuildLevelAndIndex();

    		std::cout << "New Grid Size: " << grid_->getSize() << std::endl;
    		alpha_->resizeZero(grid_->getSize());
    	}
    	else
    	{
    		std::cout << "Grid Size: " << grid_->getSize() << std::endl;
    	}

    	sg::base::DataVector b(alpha_->getSize());
    	SLESystem.generateb(classes, b);

    	if (i == AdaptConfig.noRefinements_)
    	{
    		myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
    		myCG->setEpsilon(SolverConfigFinal.eps_);
    	}
    	myCG->solve(SLESystem, *alpha_, b, true, false, 0.0);

        execTime += myStopwatch->stop();

    	std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
    	std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;

    	// Calc flops and mem bandwidth
    	//nGridsize = grid_->getSize();

    	//calcGFlopsAndGBytes(gridtype, myGrid, nInstancesNo, nGridsize, nDim, myCG->getNumberIterations(), sizeof(double), GFlops, GBytes);

    	//std::cout << std::endl;
        //std::cout << "Current GFlop/s: " << GFlops/execTime << std::endl;
        //std::cout << "Current GByte/s: " << GBytes/execTime << std::endl;
        //std::cout << std::endl;

//#ifndef TEST_LAST_ONLY
//		// Do tests on test data
//    	if (isRegression)
//    	{
//    		sg::datadriven::OperationTest* myTest = sg::op_factory::createOperationTest(*myGrid);
//			acc = myTest->testMSE(alpha, data, classes);
//			std::cout << "MSE (train): " << acc << std::endl;
//			accTest = myTest->testMSE(alpha, testData, testclasses);
//			std::cout << "MSE (test): " << accTest << std::endl << std::endl;
//			delete myTest;
//
//			if (((i > 0) && (oldAcc <= accTest)) || accTest == 0.0)
//			{
//				std::cout << "The grid is becoming worse --> stop learning" << std::endl;
//				break;
//			}
//
//			oldAcc = accTest;
//    	}
//    	else
//    	{
//    		sg::datadriven::OperationTest* myTest = sg::op_factory::createOperationTest(*myGrid);
//    		sg::base::DataVector charNumbers(4);
//    		acc = myTest->testWithCharacteristicNumber(alpha, data, classes, charNumbers);
//    		acc /= static_cast<double>(classes.getSize());
//    		std::cout << "train accuracy: " << acc << std::endl;
//    		std::cout << "train sensitivity: " << charNumbers[0]/(charNumbers[0] + charNumbers[3]) << std::endl;
//    		std::cout << "train specificity: " << charNumbers[1]/(charNumbers[1] + charNumbers[2]) << std::endl;
//    		std::cout << "train precision: " << charNumbers[0]/(charNumbers[0] + charNumbers[2]) << std::endl << std::endl;
//    		std::cout << "train true positives: " << charNumbers[0] << std::endl;
//    		std::cout << "train true negatives: " << charNumbers[1] << std::endl;
//    		std::cout << "train false positives: " << charNumbers[2] << std::endl;
//    		std::cout << "train false negatives: " << charNumbers[3] << std::endl << std::endl;
//    		accTest = myTest->testWithCharacteristicNumber(alpha, testData, testclasses, charNumbers);
//    		accTest /= static_cast<double>(testclasses.getSize());
//    		std::cout << "test accuracy: " << accTest << std::endl;
//    		std::cout << "test sensitivity: " << charNumbers[0]/(charNumbers[0] + charNumbers[3]) << std::endl;
//    		std::cout << "test specificity: " << charNumbers[1]/(charNumbers[1] + charNumbers[2]) << std::endl;
//    		std::cout << "test precision: " << charNumbers[0]/(charNumbers[0] + charNumbers[2]) << std::endl << std::endl;
//    		std::cout << "test true positives: " << charNumbers[0] << std::endl;
//    		std::cout << "test true negatives: " << charNumbers[1] << std::endl;
//    		std::cout << "test false positives: " << charNumbers[2] << std::endl;
//    		std::cout << "test false negatives: " << charNumbers[3] << std::endl << std::endl;
//#ifdef STORE_ROC_CURVE
//    		std::cout << "calculating ROC curves ..." << std::endl;
//    		sg::base::DataVector thresholds(ROC_POINTS+1);
//    		double gap = 2.0/((double)ROC_POINTS);
//    		for (size_t t = 0; t < ROC_POINTS+1; t++)
//    		{
//    			thresholds.set(t, 1.0-(t*gap));
//    		}
//    		sg::base::DataMatrix ROC_train(ROC_POINTS+1, 2);
//    		sg::base::DataMatrix ROC_test(ROC_POINTS+1, 2);
//    		myTest->calculateROCcurve(alpha, data, classes, thresholds, ROC_train);
//    		myTest->calculateROCcurve(alpha, testData, testclasses, thresholds, ROC_test);
//    		std::cout << "calculating ROC curves done!" << std::endl << std::endl;
//    		std::stringstream filetrain;
//    		filetrain << tfileTrain << "_SP" << "_level_" << start_level << "_refines_" << refine_count << "_lambda_" << lambda << ".roc";
//    		std::stringstream filetest;
//    		filetest << tfileTest << "_SP" << "_level_" << start_level << "_refines_" << refine_count << "_lambda_" << lambda << ".roc";
//        	storeROCcurve(ROC_train, filetrain.str());
//        	storeROCcurve(ROC_test, filetest.str());
//#endif
//    		delete myTest;
//
//			if (((i > 0) && (oldAcc >= accTest)) || accTest == 1.0)
//			{
//				std::cout << "The grid is becoming worse --> stop learning" << std::endl;
//				break;
//			}
//
//			oldAcc = accTest;
//    	}
//#endif
    }

    delete myStopwatch;

    std::cout << "Finished Learning!" << std::endl;

//#ifdef TEST_LAST_ONLY
//    std::cout << std::endl << std::endl;
//    if (isRegression)
//	{
//		sg::datadriven::OperationTest* myTest = sg::op_factory::createOperationTest(*myGrid);
//		acc = myTest->testMSE(alpha, data, classes);
//		std::cout << "MSE (train): " << acc << std::endl;
//		accTest = myTest->testMSE(alpha, testData, testclasses);
//		std::cout << "MSE (test): " << accTest << std::endl << std::endl;
//		delete myTest;
//	}
//	else
//	{
//		sg::datadriven::OperationTest* myTest = sg::op_factory::createOperationTest(*myGrid);
//		sg::base::DataVector charNumbers(4);
//		acc = myTest->testWithCharacteristicNumber(alpha, data, classes, charNumbers);
//		acc /= static_cast<double>(classes.getSize());
//		std::cout << "train accuracy: " << acc << std::endl;
//		std::cout << "train sensitivity: " << charNumbers[0]/(charNumbers[0] + charNumbers[3]) << std::endl;
//		std::cout << "train specificity: " << charNumbers[1]/(charNumbers[1] + charNumbers[2]) << std::endl;
//		std::cout << "train precision: " << charNumbers[0]/(charNumbers[0] + charNumbers[2]) << std::endl << std::endl;
//		std::cout << "train true positives: " << charNumbers[0] << std::endl;
//		std::cout << "train true negatives: " << charNumbers[1] << std::endl;
//		std::cout << "train false positives: " << charNumbers[2] << std::endl;
//		std::cout << "train false negatives: " << charNumbers[3] << std::endl << std::endl;
//		accTest = myTest->testWithCharacteristicNumber(alpha, testData, testclasses, charNumbers);
//		accTest /= static_cast<double>(testclasses.getSize());
//		std::cout << "test accuracy: " << accTest << std::endl;
//		std::cout << "test sensitivity: " << charNumbers[0]/(charNumbers[0] + charNumbers[3]) << std::endl;
//		std::cout << "test specificity: " << charNumbers[1]/(charNumbers[1] + charNumbers[2]) << std::endl;
//		std::cout << "test precision: " << charNumbers[0]/(charNumbers[0] + charNumbers[2]) << std::endl << std::endl;
//		std::cout << "test true positives: " << charNumbers[0] << std::endl;
//		std::cout << "test true negatives: " << charNumbers[1] << std::endl;
//		std::cout << "test false positives: " << charNumbers[2] << std::endl;
//		std::cout << "test false negatives: " << charNumbers[3] << std::endl << std::endl;
//#ifdef STORE_ROC_CURVE
//		std::cout << "calculating ROC curves ..." << std::endl;
//		sg::base::DataVector thresholds(ROC_POINTS+1);
//		double gap = 2.0/((double)ROC_POINTS);
//		for (size_t t = 0; t < ROC_POINTS+1; t++)
//		{
//			thresholds.set(t, 1.0-(t*gap));
//		}
//		sg::base::DataMatrix ROC_train(ROC_POINTS+1, 2);
//		sg::base::DataMatrix ROC_test(ROC_POINTS+1, 2);
//		myTest->calculateROCcurve(alpha, data, classes, thresholds, ROC_train);
//		myTest->calculateROCcurve(alpha, testData, testclasses, thresholds, ROC_test);
//		std::cout << "calculating ROC curves done!" << std::endl << std::endl;
//		std::stringstream filetrain;
//		filetrain << tfileTrain << "_SP" << "_level_" << start_level << "_refines_" << refine_count << "_lambda_" << lambda << ".roc";
//		std::stringstream filetest;
//		filetest << tfileTest << "_SP" << "_level_" << start_level << "_refines_" << refine_count << "_lambda_" << lambda << ".roc";
//    	storeROCcurve(ROC_train, filetrain.str());
//    	storeROCcurve(ROC_test, filetest.str());
//#endif
//    	delete myTest;
//	}
//#endif
//
//#ifdef GNUPLOT
//	if (nDim <= 2)
//	{
//		sg::base::GridPrinter* myPrinter = new sg::base::GridPrinter(*myGrid);
//		myPrinter->printGrid(alpha, "ClassifyBenchmark.gnuplot", GRDIRESOLUTION);
//		delete myPrinter;
//	}
//#endif
//
//    std::cout << std::endl;
//    std::cout << "===============================================================" << std::endl;
//    printSettings(dataFile, testFile, isRegression, start_level,
//			lambda, cg_max, cg_eps, refine_count, refine_thresh, refine_points, gridtype, vectorization, myGrid->getSize(), acc, accTest, execTime);
//#ifdef ITERATIVE
//    std::cout << "Needed time: " << execTime << " seconds (Double Precision)" << std::endl;
//    std::cout << std::endl << "Timing Details:" << std::endl;
//    double computeMult, completeMult, computeMultTrans, completeMultTrans;
//    mySystem->getTimers(completeMult, computeMult, completeMultTrans, computeMultTrans);
//    std::cout << "         mult (complete): " << completeMult << " seconds" << std::endl;
//    std::cout << "         mult (compute) : " << computeMult << " seconds" << std::endl;
//    std::cout << "  mult trans. (complete): " << completeMultTrans << " seconds" << std::endl;
//    std::cout << "  mult trans. (compute) : " << computeMultTrans << " seconds" << std::endl;
//    std::cout << std::endl << std::endl;
//    std::cout << "GFlop/s: " << GFlops/execTime << std::endl;
//    std::cout << "GByte/s: " << GBytes/execTime << std::endl;
//#else
//    std::cout << "Needed time: " << execTime << " seconds (Double Precision, recursive)" << std::endl;
//#endif
//    std::cout << "===============================================================" << std::endl;
//    std::cout << std::endl;
//
//#ifndef ITERATIVE
//    delete myC;
//#endif
//    delete myCG;
//    delete mySystem;
//    delete myGrid;
}

void Learner::trainGrid(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
		sg::solver::SLESolverConfiguration& SolverConfig, sg::base::OperationMatrix& SLESystem)
{

}

void Learner::train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
		sg::base::RegularGridConfiguration& GridConfig, sg::base::AdpativityConfiguration& AdaptConfig,
		sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal)
{

}

void Learner::train(sg::base::DataMatrix& testDataset, sg::base::RegularGridConfiguration& GridConfig,
		sg::solver::SLESolverConfiguration& SolverConfig)
{

}

void Learner::test(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes, bool isRegression)
{

}

void Learner::store(std::string tGridFilename, std::string tAlphaFilename)
{

}

sg::base::DataVector Learner::getAlpha()
{
	sg::base::DataVector a(1);
	return a;
}


}

}
