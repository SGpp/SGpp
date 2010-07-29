/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBLinear.hpp"
#include "tools/datadriven/ARFFTools.hpp"
#include "common/AlignedMemory.hpp"

#include <string>
#include <iostream>

#define LEVELS 6
#define DATAFILE "DR5_nowarnings_less05_train.arff"
#define ITERATIONS 1

/**
 * executes the OperationB mult_transposed (multiple function evaluations)
 * with the DR5 testcase
 */
void executesOperationBmultTrans_DR5()
{
    // Write the data of CG
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Performance Test App" << std::endl;

	double execTime = 0.0;
	sg::ARFFTools ARFFTool;
	std::string tfileTrain = DATAFILE;

	size_t nDim = ARFFTool.getDimension(tfileTrain);
	size_t nInstancesNo = ARFFTool.getNumberInstances(tfileTrain);

	// Create Grid
	sg::Grid* myGrid = new sg::LinearGrid(nDim);

	// Generate regular Grid with LEVELS Levels
	sg::GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(LEVELS);
	delete myGenerator;

	// Read data from file
	DataVector data(nInstancesNo, nDim);
    DataVector result(nInstancesNo);
    DataVector alpha(myGrid->getSize());

    // Set DataVectors
    ARFFTool.readTrainingData(tfileTrain, data);
    result.setAll(0.0);
    alpha.setAll(1.0);

    // init the Systemmatrix Functor
    sg::OperationBLinear* B = (sg::OperationBLinear*)myGrid->createOperationB();

    std::cout << "Grid, Data, Operation created! Start test..." << std::endl;
    std::cout << "GridSize: " << myGrid->getSize() << std::endl;
    std::cout << "DataSize: " << data.getSize() << std::endl;

    // Generate SOA from AOS
    double* level = new double[myGrid->getSize()*nDim];
    double* index = new double[myGrid->getSize()*nDim];

    sg::SGppStopwatch* myStopwatch = new sg::SGppStopwatch();
    myStopwatch->start();

    myGrid->getStorage()->getLevelIndexArraysForEval(level, index);

    for (size_t i = 0; i < ITERATIONS; i++)
    {
    	//B->multTranspose(alpha, data, result);
    	B->multTransposeIterative(level,index,alpha, data, result);
    }

    execTime = myStopwatch->stop();

    delete[] level;
    delete[] index;
    delete B;

    // Write the data of CG
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Needed time: " << execTime << " seconds" << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << std::endl;

    delete myGrid;
}


/**
 * Testapplication for the Intel VTune Profiling Tool
 * and a measurement app for Sparse Grid Algorithms building blocks
 */
int main(int argc, char *argv[])
{
	executesOperationBmultTrans_DR5();

	return 0;
}
