/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "tools/datadriven/ARFFTools.hpp"

#include <string>
#include <iostream>

#define DATAFILE "DR5_nowarnings_less05_train.arff"
//#define DATAFILE "twospirals.wieland.arff"

//#define TESTFILE "DR5_nowarnings_less05_test.arff"

#define LEVELS 3
#define ITERATIONS 10
#define REFINEMENTS 1
#define CG_IMAX 1000
#define CG_EPS 0.000001
#define LAMBDA 0.00001
#define REFINE_THRESHOLD 0.0
#define REFINE_NUM 100

#define USE_SSE


void adaptRegressionTest()
{
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Regression/Classification Test App" << std::endl;
    std::cout << "===============================================================" << std::endl;

	double execTime = 0.0;
	sg::ARFFTools ARFFTool;
	std::string tfileTrain = DATAFILE;
//	std::string tfileTest = TESTFILE;

	size_t nDim = ARFFTool.getDimension(tfileTrain);
	size_t nInstancesNo = ARFFTool.getNumberInstances(tfileTrain);

	// Create Grid
	sg::Grid* myGrid = new sg::LinearGrid(nDim);

	// Generate regular Grid with LEVELS Levels
	sg::GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(LEVELS);
	delete myGenerator;

	// Read data from file
	DataMatrix data(nInstancesNo, nDim);
    DataVector classes(nInstancesNo);
//    DataVector testclasses(nInstancesNo);
    DataVector result(nInstancesNo);
    DataVector alpha(myGrid->getSize());

    // Set DataVectors
    ARFFTool.readTrainingData(tfileTrain, data);
    ARFFTool.readClasses(tfileTrain, classes);
//    ARFFTool.readClasses(tfileTest, testclasses);
    result.setAll(0.0);
    alpha.setAll(0.0);

    // Generate CG to solve System
    sg::ConjugateGradients* myCG = new sg::ConjugateGradients(CG_IMAX, CG_EPS);
#ifdef USE_SSE
    sg::DMSystemMatrixSSEIdentity* mySystem = new sg::DMSystemMatrixSSEIdentity(*myGrid, data, LAMBDA);
#else
    sg::OperationMatrix* myC = myGrid->createOperationIdentity();
    sg::DMSystemMatrix* mySystem = new sg::DMSystemMatrix(*myGrid, data, *myC, LAMBDA);
#endif

    std::cout << "Starting Learning...." << std::endl;
    // execute adaptsteps
    sg::SGppStopwatch* myStopwatch = new sg::SGppStopwatch();
    myStopwatch->start();
    for (size_t i = 0; i < REFINEMENTS; i++)
    {
    	std::cout << "Doing refinement :" << i << std::endl;

    	// Do Refinements
    	if (i > 0)
    	{
    		sg::SurplusRefinementFunctor* myRefineFunc = new sg::SurplusRefinementFunctor(&alpha, REFINE_NUM, REFINE_THRESHOLD);
    		myGrid->createGridGenerator()->refine(myRefineFunc);
    		delete myRefineFunc;
#ifdef USE_SSE
    		mySystem->rebuildLevelAndIndex();
#endif
    		std::cout << "New Grid Size: " << myGrid->getStorage()->size() << std::endl;
    		alpha.resizeZero(myGrid->getStorage()->size());
    	}

    	DataVector b(alpha.getSize());
    	mySystem->generateb(classes, b);

    	//std::cout << b.toString() << std::endl << std::endl;
    	//std::cout << alpha.toString() << std::endl << std::endl;

    	myCG->solve(*mySystem, alpha, b, false, false, 0.0);

    	std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
    	std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;

    	// Do tests on test data
//        sg::OperationTest* myTest = myGrid->createOperationTest();
//        correct = myTest->test(alpha, data, testclasses);
//        delete myTest;
    }
    execTime = myStopwatch->stop();
    delete myStopwatch;

    std::cout << "Finished Learning!" << std::endl;

//    sg::GridPrinter* myPrinter = new sg::GridPrinter(*myGrid);
//    myPrinter->printGrid(alpha, "VTuneTest_Result.gnuplot", 50);
//    delete myPrinter;

    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Needed time: " << execTime << " seconds" << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << std::endl;

#ifndef USE_SSE
    delete myC;
#endif
    delete myCG;
    delete mySystem;
    delete myGrid;
}

/**
 * Testapplication for the Intel VTune Profiling Tool
 * and a measurement app for Sparse Grid Algorithms building blocks
 */
int main(int argc, char *argv[])
{
	adaptRegressionTest();

	return 0;
}
