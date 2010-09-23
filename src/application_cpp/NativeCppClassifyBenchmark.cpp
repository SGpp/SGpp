/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "data/DataVectorSP.hpp"
#include "data/DataMatrixSP.hpp"
#include "algorithm/datadriven/DMSystemMatrixSPVectorizedIdentity.hpp"
#include "solver/sle/ConjugateGradientsSP.hpp"
#include "tools/datadriven/ARFFTools.hpp"

#include <string>
#include <iostream>

//#define DATAFILE "DR5_nowarnings_less05_train.arff"
//#define DATAFILE "twospirals.wieland.arff"
//#define DATAFILE "liver-disorders_normalized.arff"
//#define DATAFILE "ripleyGarcke.train.arff"
#define DATAFILE "chess_02D_tr.dat.arff"
//#define DATAFILE "chess_05D_tr.dat.arff"

//#define TESTFILE "DR5_nowarnings_less05_test.arff"
//#define TESTFILE "twospirals.wieland.arff"
//#define TESTFILE "liver-disorders_normalized.arff"
//#define TESTFILE "ripleyGarcke.test.arff"
#define TESTFILE "chess_02D_te.dat.arff"
//#define TESTFILE "chess_05D_te.dat.arff"

// grid generation settingsd
#define LEVELS 8
#define REFINEMENTS 0
#define REFINE_THRESHOLD 0.0
#define REFINE_NUM_POINTS 100

// solving settings
#define CG_IMAX 10000
#define CG_EPS 0.0001

// regularization fector
#define LAMBDA 0.0001

// print grid in gnuplot readable format (1D and 2D only)
#define GNUPLOT
#define GRDIRESOLUTION 50

// at least one has to be defined, otherwise scalar&recursive version is used for DP, SSE for SP
//#define USE_SSE
//#define USE_AVX
#define USE_OCL

// define if you want to use single precision floats (may deliver speed-up of 2 or greater),
// BUT: CG method may not converge because of bad system matrix condition.
#define USEFLOAT

// define this if you want to execute a regression
//#define EXEC_REGRESSION

// define this if you want to use grids with Neumann boundaries.
#define USE_BOUNDARIES

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

void printSettings()
{
	std::cout << std::endl;
	std::cout << "Train dataset: " << DATAFILE << std::endl;
	std::cout << "Test dataset: " << TESTFILE << std::endl;
	std::cout << "Startlevel: " << LEVELS << std::endl << std::endl;

	std::cout << "Num. Refinements: " << REFINEMENTS << std::endl;
	std::cout << "Refine Threshold: " << REFINE_THRESHOLD << std::endl;
	std::cout << "Refine number points: " << REFINE_NUM_POINTS << std::endl << std::endl;

	std::cout << "Max. CG Iterations: " << CG_IMAX << std::endl;
	std::cout << "CG epsilon: " << CG_EPS << std::endl << std::endl;

	std::cout << "Lambda: " << LAMBDA << std::endl;

#ifdef USE_SSE
	std::cout << "Vectorized: SSE" << std::endl << std::endl;
#endif
#ifdef USE_AVX
	std::cout << "Vectorized: AVX" << std::endl << std::endl;
#endif
#ifdef USE_OCL
	std::cout << "Vectorized: OpenCL (nVidia Fermi tested)" << std::endl << std::endl;
#endif

#ifdef EXEC_REGRESSION
	std::cout << "Mode: Regression" << std::endl << std::endl;
#else
	std::cout << "Mode: Classification" << std::endl << std::endl;
#endif

#ifdef USE_BOUNDARIES
	std::cout << "Boundary-Mode: Neumann" << std::endl << std::endl;
#else
	std::cout << "Boundary-Mode: Dirichlet-0" << std::endl << std::endl;
#endif
}

void adaptClassificationTest(bool isRegression)
{
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_OCL)
    std::cout << "Classification Test App (Double Precision)" << std::endl;
#else
    std::cout << "Classification Test App (Double Precision, recursive)" << std::endl;
#endif
    std::cout << "===============================================================" << std::endl << std::endl;

    printSettings();

    double execTime = 0.0;
	sg::ARFFTools ARFFTool;
	std::string tfileTrain = DATAFILE;
	std::string tfileTest = TESTFILE;

	size_t nDim = ARFFTool.getDimension(tfileTrain);
	size_t nInstancesNo = ARFFTool.getNumberInstances(tfileTrain);
	size_t nInstancesTestNo = ARFFTool.getNumberInstances(tfileTest);

	// Create Grid
	sg::Grid* myGrid;
#ifdef USE_BOUNDARIES
	myGrid = new sg::LinearTrapezoidBoundaryGrid(nDim);
#else
	myGrid = new sg::LinearGrid(nDim);
#endif

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

    // Variable to save MSE/Acc from former iteration
    double oldAcc = 0.0;

    // Generate CG to solve System
    sg::ConjugateGradients* myCG = new sg::ConjugateGradients(CG_IMAX, CG_EPS);
#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_OCL)
#ifdef USE_SSE
    sg::DMSystemMatrixVectorizedIdentity* mySystem = new sg::DMSystemMatrixVectorizedIdentity(*myGrid, data, LAMBDA, "SSE");
#endif
#ifdef USE_AVX
    sg::DMSystemMatrixVectorizedIdentity* mySystem = new sg::DMSystemMatrixVectorizedIdentity(*myGrid, data, LAMBDA, "AVX");
#endif
#ifdef USE_OCL
    sg::DMSystemMatrixVectorizedIdentity* mySystem = new sg::DMSystemMatrixVectorizedIdentity(*myGrid, data, LAMBDA, "OCL");
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
    	std::cout << std::endl << "Doing refinement :" << i << std::endl;

    	// Do Refinements
    	if (i > 0)
    	{
    		sg::SurplusRefinementFunctor* myRefineFunc = new sg::SurplusRefinementFunctor(&alpha, REFINE_NUM_POINTS, REFINE_THRESHOLD);
    		myGrid->createGridGenerator()->refine(myRefineFunc);
    		delete myRefineFunc;

#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_OCL)
    		mySystem->rebuildLevelAndIndex();
#endif

    		std::cout << "New Grid Size: " << myGrid->getStorage()->size() << std::endl;
    		alpha.resizeZero(myGrid->getStorage()->size());
    	}
    	else
    	{
    		std::cout << "Grid Size: " << myGrid->getStorage()->size() << std::endl;
    	}

    	DataVector b(alpha.getSize());
    	mySystem->generateb(classes, b);

    	myCG->solve(*mySystem, alpha, b, true, true, 0.0);

    	std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
    	std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;

		// Do tests on test data
    	if (isRegression)
    	{
    		sg::OperationTest* myTest = myGrid->createOperationTest();
			double mse = myTest->testMSE(alpha, data, classes);
			std::cout << "MSE (train): " << mse << std::endl;
			double mseTest = myTest->testMSE(alpha, testData, testclasses);
			std::cout << "MSE (test): " << mseTest << std::endl;
			delete myTest;

			if (((i > 0) && (oldAcc <= mseTest)) || mseTest == 0.0)
			{
				std::cout << "The grid is becoming worse --> stop learning" << std::endl;
				break;
			}

			oldAcc = mseTest;
    	}
    	else
    	{
    		sg::OperationTest* myTest = myGrid->createOperationTest();
			double acc = myTest->test(alpha, data, classes);
			acc /= static_cast<double>(classes.getSize());
			std::cout << "train acc.: " << acc << std::endl;
			double accTest = myTest->test(alpha, testData, testclasses);
			accTest /= static_cast<double>(testclasses.getSize());
			std::cout << "test acc.: " << accTest << std::endl;
			delete myTest;

			if (((i > 0) && (oldAcc >= accTest)) || accTest == 1.0)
			{
				std::cout << "The grid is becoming worse --> stop learning" << std::endl;
				break;
			}

			oldAcc = accTest;
    	}
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
    printSettings();
#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_OCL)
    std::cout << "Needed time: " << execTime << " seconds (Double Precision)" << std::endl;
    std::cout << std::endl << "Timing Details:" << std::endl;
    double computeMult, completeMult, computeMultTrans, completeMultTrans;
    mySystem->getTimers(completeMult, computeMult, completeMultTrans, computeMultTrans);
    std::cout << "         mult (complete): " << completeMult << " seconds" << std::endl;
    std::cout << "         mult (compute) : " << computeMult << " seconds" << std::endl;
    std::cout << "  mult trans. (complete): " << completeMultTrans << " seconds" << std::endl;
    std::cout << "  mult trans. (compute) : " << computeMultTrans << " seconds" << std::endl;
#else
    std::cout << "Needed time: " << execTime << " seconds (Double Precision, recursive)" << std::endl;
#endif
    std::cout << "===============================================================" << std::endl;
    std::cout << std::endl;

#ifndef USE_SSE
#ifndef USE_AVX
#ifndef USE_OCL
    delete myC;
#endif
#endif
#endif
    delete myCG;
    delete mySystem;
    delete myGrid;
}


void adaptClassificationTestSP(bool isRegression)
{
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Classification Test App (Single Precision)" << std::endl;
    std::cout << "===============================================================" << std::endl << std::endl;

    printSettings();

	double execTime = 0.0;
	sg::ARFFTools ARFFTool;
	std::string tfileTrain = DATAFILE;
	std::string tfileTest = TESTFILE;

	size_t nDim = ARFFTool.getDimension(tfileTrain);
	size_t nInstancesNo = ARFFTool.getNumberInstances(tfileTrain);
	size_t nInstancesTestNo = ARFFTool.getNumberInstances(tfileTest);

	// Create Grid
	sg::Grid* myGrid;
#ifdef USE_BOUNDARIES
	myGrid = new sg::LinearTrapezoidBoundaryGrid(nDim);
#else
	myGrid = new sg::LinearGrid(nDim);
#endif

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

    // Variable to save MSE/Acc from former iteration
    double oldAcc = 0.0;

    // Generate CG to solve System
    sg::ConjugateGradientsSP* myCG = new sg::ConjugateGradientsSP(CG_IMAX, CG_EPS);

#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_OCL)
#ifdef USE_SSE
    sg::DMSystemMatrixSPVectorizedIdentity* mySystem = new sg::DMSystemMatrixSPVectorizedIdentity(*myGrid, dataSP, LAMBDA, "SSE");
#endif
#ifdef USE_AVX
    sg::DMSystemMatrixSPVectorizedIdentity* mySystem = new sg::DMSystemMatrixSPVectorizedIdentity(*myGrid, dataSP, LAMBDA, "AVX");
#endif
#ifdef USE_OCL
    sg::DMSystemMatrixSPVectorizedIdentity* mySystem = new sg::DMSystemMatrixSPVectorizedIdentity(*myGrid, dataSP, LAMBDA, "OCL");
#endif
#else
    sg::DMSystemMatrixSPVectorizedIdentity* mySystem = new sg::DMSystemMatrixSPVectorizedIdentity(*myGrid, dataSP, LAMBDA, "SSE");
#endif

    std::cout << "Starting Learning...." << std::endl;
    // execute adaptsteps
    sg::SGppStopwatch* myStopwatch = new sg::SGppStopwatch();
    myStopwatch->start();
    for (size_t i = 0; i < REFINEMENTS+1; i++)
    {
    	std::cout << std::endl << "Doing refinement :" << i << std::endl;

    	// Do Refinements
    	if (i > 0)
    	{
    		convertDataVectorSPToDataVector(alphaSP, alpha);
    		sg::SurplusRefinementFunctor* myRefineFunc = new sg::SurplusRefinementFunctor(&alpha, REFINE_NUM_POINTS, REFINE_THRESHOLD);
    		myGrid->createGridGenerator()->refine(myRefineFunc);
    		delete myRefineFunc;

#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_OCL)
    		mySystem->rebuildLevelAndIndex();
#endif

    		std::cout << "New Grid Size: " << myGrid->getStorage()->size() << std::endl;
    		alpha.resizeZero(myGrid->getStorage()->size());
    		alphaSP.resizeZero(myGrid->getStorage()->size());
    	}
    	else
    	{
    		std::cout << "Grid Size: " << myGrid->getStorage()->size() << std::endl;
    	}

    	DataVectorSP bSP(alphaSP.getSize());
    	mySystem->generateb(classesSP, bSP);

    	myCG->solve(*mySystem, alphaSP, bSP, true, true, 0.0);

    	std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
    	std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;

    	// Do tests on test data
    	convertDataVectorSPToDataVector(alphaSP, alpha);
    	if (isRegression)
    	{
    		sg::OperationTest* myTest = myGrid->createOperationTest();
			double mse = myTest->testMSE(alpha, data, classes);
			std::cout << "MSE (train): " << mse << std::endl;
			double mseTest = myTest->testMSE(alpha, testData, testclasses);
			std::cout << "MSE (test): " << mseTest << std::endl;
			delete myTest;

			if (((i > 0) && (oldAcc <= mseTest)) || mseTest == 0.0)
			{
				std::cout << "The grid is becoming worse --> stop learning" << std::endl;
				break;
			}

			oldAcc = mseTest;
    	}
    	else
    	{
    		sg::OperationTest* myTest = myGrid->createOperationTest();
			double acc = myTest->test(alpha, data, classes);
			acc /= static_cast<double>(classes.getSize());
			std::cout << "train acc.: " << acc << std::endl;
			double accTest = myTest->test(alpha, testData, testclasses);
			accTest /= static_cast<double>(testclasses.getSize());
			std::cout << "test acc.: " << accTest << std::endl;
			delete myTest;

			if (((i > 0) && (oldAcc >= accTest)) || accTest == 1.0)
			{
				std::cout << "The grid is becoming worse --> stop learning" << std::endl;
				break;
			}

			oldAcc = accTest;
    	}
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
    printSettings();
    std::cout << "Needed time: " << execTime << " seconds (Single Precision)" << std::endl;
    std::cout << std::endl << "Timing Details:" << std::endl;
    double computeMult, completeMult, computeMultTrans, completeMultTrans;
    mySystem->getTimers(completeMult, computeMult, completeMultTrans, computeMultTrans);
    std::cout << "         mult (complete): " << completeMult << " seconds" << std::endl;
    std::cout << "         mult (compute) : " << computeMult << " seconds" << std::endl;
    std::cout << "  mult trans. (complete): " << completeMultTrans << " seconds" << std::endl;
    std::cout << "  mult trans. (compute) : " << computeMultTrans << " seconds" << std::endl;
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
#ifdef EXEC_REGRESSION
	adaptClassificationTestSP(true);
#else
	adaptClassificationTestSP(false);
#endif
#else
#ifdef EXEC_REGRESSION
	adaptClassificationTest(true);
#else
	adaptClassificationTest(false);
#endif
#endif
	return 0;
}
