/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "data/DataVectorSP.hpp"
#include "data/DataMatrixSP.hpp"
#include "algorithm/datadriven/DMSystemMatrixSPSSEIdentity.hpp"
#include "solver/sle/ConjugateGradientsSP.hpp"
#include "tools/datadriven/ARFFTools.hpp"

#include <string>
#include <iostream>

//#define DATAFILE "DR5_nowarnings_less05_train.arff"
//#define DATAFILE "twospirals.wieland.arff"
//#define DATAFILE "liver-disorders_normalized.arff"
#define DATAFILE "ripleyGarcke.train.arff"

//#define TESTFILE "DR5_nowarnings_less05_test.arff"
//#define TESTFILE "twospirals.wieland.arff"
//#define TESTFILE "liver-disorders_normalized.arff"
#define TESTFILE "ripleyGarcke.test.arff"

// grid generation settings
#define LEVELS 3
#define REFINEMENTS 2
#define REFINE_THRESHOLD 0.0
#define REFINE_NUM_POINTS 3

// solving settings
#define CG_IMAX 10000
#define CG_EPS 0.00001

// regularization fector
#define LAMBDA 0.000001

// tests learned grid with test data
#define TESTRESULT

// print grid in gnuplot readable format (1D and 2D only)
#define GNUPLOT
#define GRDIRESOLUTION 50

// at least one has to be defined
#define USE_SSE
//#define USE_AVX

// define if you want to use single precision floats (may deliver speed-up of 2 or great),
// BUT: CG method may not converge because of bad system matrix condition.
//#define USEFLOAT

void convertDataVectorToDataVectorSP(DataVector& src, DataVectorSP& dest)
{
	if (src.getSize() != dest.getSize())
	{
		return;
	}
	else
	{
		for (size_t i = 0; i < src.getSize(); i++)
		{
			dest.set(i, static_cast<float>(src.get(i)));
		}
	}
}

void convertDataMatrixToDataMatrixSP(DataMatrix& src, DataMatrixSP& dest)
{
	if (src.getNcols() != dest.getNcols() || src.getNrows() != dest.getNrows())
	{
		return;
	}
	else
	{
		for (size_t i = 0; i < src.getNrows(); i++)
		{
			for (size_t j = 0; j < src.getNcols(); j++)
			{
				dest.set(i, j, static_cast<float>(src.get(i, j)));
			}
		}
	}
}

void convertDataVectorSPToDataVector(DataVectorSP& src, DataVector& dest)
{
	if (src.getSize() != dest.getSize())
	{
		return;
	}
	else
	{
		for (size_t i = 0; i < src.getSize(); i++)
		{
			dest.set(i, static_cast<double>(src.get(i)));
		}
	}
}

void adaptRegressionTest()
{
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Regression/Classification Test App (Double Precision)" << std::endl;
    std::cout << "===============================================================" << std::endl << std::endl;

	double execTime = 0.0;
	sg::ARFFTools ARFFTool;
	std::string tfileTrain = DATAFILE;
	std::string tfileTest = TESTFILE;

	size_t nDim = ARFFTool.getDimension(tfileTrain);
	size_t nInstancesNo = ARFFTool.getNumberInstances(tfileTrain);
	size_t nInstancesTestNo = ARFFTool.getNumberInstances(tfileTest);

	// Create Grid
	sg::Grid* myGrid;

	myGrid = new sg::LinearGrid(nDim);
	// Generate regular Grid with LEVELS Levels
	sg::GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(LEVELS);
	delete myGenerator;

	// Define DP data
	DataMatrix data(nInstancesNo, nDim);
    DataVector classes(nInstancesNo);
    DataMatrix testData(nInstancesTestNo, nDim);
    DataVector testclasses(nInstancesTestNo);
    DataVector result(nInstancesNo);
    DataVector alpha(myGrid->getSize());

	// Read data from file
    ARFFTool.readTrainingData(tfileTrain, data);
    ARFFTool.readTrainingData(tfileTest, testData);
    ARFFTool.readClasses(tfileTrain, classes);
    ARFFTool.readClasses(tfileTest, testclasses);

    result.setAll(0.0);
    alpha.setAll(0.0);

    // Generate CG to solve System
    sg::ConjugateGradients* myCG = new sg::ConjugateGradients(CG_IMAX, CG_EPS);
#if defined(USE_SSE) || defined(USE_AVX)
#ifdef USE_SSE
    sg::DMSystemMatrixSSEIdentity* mySystem = new sg::DMSystemMatrixSSEIdentity(*myGrid, data, LAMBDA);
#endif
#ifdef USE_AVX
    sg::DMSystemMatrixAVXIdentity* mySystem = new sg::DMSystemMatrixAVXIdentity(*myGrid, data, LAMBDA);
#endif
#else
    sg::OperationMatrix* myC = myGrid->createOperationIdentity();
    //sg::OperationMatrix* myC = myGrid->createOperationLaplace();
    sg::DMSystemMatrix* mySystem = new sg::DMSystemMatrix(*myGrid, data, *myC, LAMBDA);
#endif

    std::cout << "Starting Learning...." << std::endl;
    // execute adaptsteps
    sg::SGppStopwatch* myStopwatch = new sg::SGppStopwatch();
    myStopwatch->start();
    for (size_t i = 0; i < REFINEMENTS+1; i++)
    {
    	std::cout << "Doing refinement :" << i << std::endl;

    	// Do Refinements
    	if (i > 0)
    	{
    		sg::SurplusRefinementFunctor* myRefineFunc = new sg::SurplusRefinementFunctor(&alpha, REFINE_NUM_POINTS, REFINE_THRESHOLD);
    		myGrid->createGridGenerator()->refine(myRefineFunc);
    		delete myRefineFunc;

#if defined(USE_SSE) || defined(USE_AVX)
    		mySystem->rebuildLevelAndIndex();
#endif

    		std::cout << "New Grid Size: " << myGrid->getStorage()->size() << std::endl;
    		alpha.resizeZero(myGrid->getStorage()->size());
    	}

    	DataVector b(alpha.getSize());
    	mySystem->generateb(classes, b);

    	myCG->solve(*mySystem, alpha, b, true, true, 0.0);

    	std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
    	std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;

#ifdef TESTRESULT
    	// Do tests on test data
        sg::OperationTest* myTest = myGrid->createOperationTest();
        double correct = myTest->test(alpha, testData, testclasses);
        std::cout << "Final test acc.: " << correct/static_cast<double>(testclasses.getSize()) << std::endl;
        delete myTest;
#endif
    }

    execTime = myStopwatch->stop();
    delete myStopwatch;

    std::cout << "Finished Learning!" << std::endl;
#ifdef GNUPLOT
	if (nDim <= 2)
	{
		sg::GridPrinter* myPrinter = new sg::GridPrinter(*myGrid);
		myPrinter->printGrid(alpha, "ClassifyBenchmark.gnuplot", GRDIRESOLUTION);
		delete myPrinter;
	}
#endif

    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Needed time: " << execTime << " seconds (Double Precision)" << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << std::endl;

#ifndef USE_SSE
#ifndef USE_AVX
    delete myC;
#endif
#endif
    delete myCG;
    delete mySystem;
    delete myGrid;
}


void adaptRegressionTestSP()
{
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Regression/Classification Test App (Single Precision)" << std::endl;
    std::cout << "===============================================================" << std::endl << std::endl;

	double execTime = 0.0;
	sg::ARFFTools ARFFTool;
	std::string tfileTrain = DATAFILE;
	std::string tfileTest = TESTFILE;

	size_t nDim = ARFFTool.getDimension(tfileTrain);
	size_t nInstancesNo = ARFFTool.getNumberInstances(tfileTrain);
	size_t nInstancesTestNo = ARFFTool.getNumberInstances(tfileTest);

	// Create Grid
	sg::Grid* myGrid;

	myGrid = new sg::LinearGrid(nDim);
	// Generate regular Grid with LEVELS Levels
	sg::GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(LEVELS);
	delete myGenerator;

	// Define DP data
	DataMatrix data(nInstancesNo, nDim);
    DataVector classes(nInstancesNo);
    DataMatrix testData(nInstancesTestNo, nDim);
    DataVector testclasses(nInstancesTestNo);
    DataVector result(nInstancesNo);
    DataVector alpha(myGrid->getSize());

    // Define SP data
	DataMatrixSP dataSP(nInstancesNo, nDim);
    DataVectorSP classesSP(nInstancesNo);
    DataVectorSP resultSP(nInstancesNo);
    DataVectorSP alphaSP(myGrid->getSize());

	// Read data from file
    ARFFTool.readTrainingData(tfileTrain, data);
    ARFFTool.readTrainingData(tfileTest, testData);
    ARFFTool.readClasses(tfileTrain, classes);
    ARFFTool.readClasses(tfileTest, testclasses);

    result.setAll(0.0);
    alpha.setAll(0.0);

    convertDataMatrixToDataMatrixSP(data, dataSP);
    convertDataVectorToDataVectorSP(alpha, alphaSP);
    convertDataVectorToDataVectorSP(classes, classesSP);
    convertDataVectorToDataVectorSP(result, resultSP);
    convertDataVectorToDataVectorSP(alpha, alphaSP);

    // Generate CG to solve System
    sg::ConjugateGradientsSP* myCG = new sg::ConjugateGradientsSP(CG_IMAX, CG_EPS);

#if defined(USE_SSE) || defined(USE_AVX)
#ifdef USE_SSE
    sg::DMSystemMatrixSPSSEIdentity* mySystem = new sg::DMSystemMatrixSPSSEIdentity(*myGrid, dataSP, LAMBDA);
#endif
#ifdef USE_AVX
    //sg::DMSystemMatrixSPAVXIdentity* mySystem = new sg::DMSystemMatrixSPAVXIdentity(*myGrid, data, LAMBDA);
#endif
#else
    sg::DMSystemMatrixSPSSEIdentity* mySystem = new sg::DMSystemMatrixSPSSEIdentity(*myGrid, dataSP, LAMBDA);
#endif

    std::cout << "Starting Learning...." << std::endl;
    // execute adaptsteps
    sg::SGppStopwatch* myStopwatch = new sg::SGppStopwatch();
    myStopwatch->start();
    for (size_t i = 0; i < REFINEMENTS+1; i++)
    {
    	std::cout << "Doing refinement :" << i << std::endl;

    	// Do Refinements
    	if (i > 0)
    	{
    		convertDataVectorSPToDataVector(alphaSP, alpha);
    		sg::SurplusRefinementFunctor* myRefineFunc = new sg::SurplusRefinementFunctor(&alpha, REFINE_NUM_POINTS, REFINE_THRESHOLD);
    		myGrid->createGridGenerator()->refine(myRefineFunc);
    		delete myRefineFunc;

#if defined(USE_SSE) || defined(USE_AVX)
    		mySystem->rebuildLevelAndIndex();
#endif

    		std::cout << "New Grid Size: " << myGrid->getStorage()->size() << std::endl;
    		alpha.resizeZero(myGrid->getStorage()->size());
    		alphaSP.resizeZero(myGrid->getStorage()->size());
    	}

    	DataVectorSP bSP(alphaSP.getSize());
    	mySystem->generateb(classesSP, bSP);

    	myCG->solve(*mySystem, alphaSP, bSP, true, true, 0.0);

    	std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
    	std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;

#ifdef TESTRESULT
    	// Do tests on test data
    	convertDataVectorSPToDataVector(alphaSP, alpha);
        sg::OperationTest* myTest = myGrid->createOperationTest();
        double correct = myTest->test(alpha, testData, testclasses);
        std::cout << "Final test acc.: " << correct/static_cast<double>(testclasses.getSize()) << std::endl;
        delete myTest;
#endif
    }

    execTime = myStopwatch->stop();
    delete myStopwatch;

    std::cout << "Finished Learning!" << std::endl;

#ifdef GNUPLOT
	if (nDim <= 2)
	{
		convertDataVectorSPToDataVector(alphaSP, alpha);
		sg::GridPrinter* myPrinter = new sg::GridPrinter(*myGrid);
		myPrinter->printGrid(alpha, "ClassifyBenchmark.gnuplot", GRDIRESOLUTION);
		delete myPrinter;
	}
#endif

    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Needed time: " << execTime << " seconds (Single Precision)" << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << std::endl;

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
#ifdef USEFLOAT
	adaptRegressionTestSP();
#else
	adaptRegressionTest();
#endif
	return 0;
}
