/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "application/classification/Classifier.hpp"
#include "exception/operation_exception.hpp"
#include "algorithm/datadriven/DMSystemMatrix.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "tools/datadriven/ARFFTools.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"
#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/type/ModLinearGrid.hpp"
#include "tools/common/SGppStopwatch.hpp"
#include <iostream>

namespace sg
{

Classifier::Classifier()
{
	this->myGrid = NULL;
	this->levels = 0;
	this->GridType = "N";
	this->dim = 0;
	this->instancesNo = 0;
	this->testinstancesNo = 0;
	this->lambda = 0.0;
	this->StiffnessMode = "L";
	this->IterationMax = 1000;
	this->epsilon = 0.0;
	std::cout << "A new Classifier instance wass created" << std::endl;
}

Classifier::~Classifier()
{
	delete myGrid;
}

void Classifier::trainNtestRegular(std::string tfileTrain, std::string tfileTest, size_t level, double lambda, std::string GridType, std::string StiffMode, double epsilon, size_t imax)
{
	std::cout << "trainNtestRegular was called within this classifier instance" << std::endl;
	ARFFTools ARFFTool;

	this->levels = level;
	this->GridType = GridType;
	this->lambda = lambda;
	this->StiffnessMode = StiffMode;
	this->IterationMax = imax;
	this->epsilon = epsilon;
	this->dim = ARFFTool.getDimension(tfileTrain);
	this->instancesNo = ARFFTool.getNumberInstances(tfileTrain);
	this->testinstancesNo = ARFFTool.getNumberInstances(tfileTest);

	std::cout << "the problem has " << this->dim << " dimensions" << std::endl;
	std::cout << "the training set contains " << this->instancesNo << " instances" << std::endl;
	std::cout << "the test set contains " << this->testinstancesNo << " instances" << std::endl;
	if (this->StiffnessMode == "L")
	{
		std::cout << "using Laplacian matrix" << std::endl;
	}
	else if (this->StiffnessMode == "I")
	{
		std::cout << "using Identity matrix" << std::endl;
	}
	else
	{
		std::cout << "using Laplacian matrix" << std::endl;
	}
	std::cout << std::endl;
	std::cout << "classifier instance is initialized by trainNtestRegular" << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "start constructing regular grid" << std::endl;

	// create the grid
	createRegularGrid();

	std::cout << "the grid has " << myGrid->getStorage()->size() << " gridpoints" << std::endl;
	std::cout << "finished construction regular grid" << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "start training grid" << std::endl;

	DataVector alpha(this->myGrid->getStorage()->size());
	alpha.setAll(0.0);

	// start the train process
	trainGrid(alpha, tfileTrain);

	std::cout << "finished training grid" << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "start testing trained grid" << std::endl;

	// start test process
	applyTestdata(alpha, tfileTest);

	std::cout << "finished testing trained grid" << std::endl;
}

void Classifier::createRegularGrid()
{
	if (this->GridType == "N")
	{
		myGrid = new LinearGrid(this->dim);
		std::cout << "A LinearGrid was created" << std::endl;
	}
	else if (this->GridType == "U")
	{
		myGrid = new LinearTrapezoidBoundaryGrid(this->dim);
		std::cout << "A LinearTrapezoidBoundaryGrid was created" << std::endl;
	}
	else if (this->GridType == "B")
	{
		myGrid = new LinearBoundaryGrid(this->dim);
		std::cout << "A LinearBoundaryGrid was created" << std::endl;
	}
	else if (this->GridType == "E")
	{
		myGrid = new ModLinearGrid(this->dim);
		std::cout << "A ModlinearGrid was created" << std::endl;
	}
	else
	{
		throw new operation_exception("No valid GridType was specified!");
	}

	GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(this->levels);
	delete myGenerator;

	std::cout << levels << " levels were added to the above created grid" << std::endl;
}

void Classifier::trainGrid(DataVector& alpha, std::string tfileTrain)
{
    ARFFTools ARFFTool;
    DataMatrix training(this->instancesNo, this->dim);
    DataVector classes(this->instancesNo);
    DataVector rhs(this->myGrid->getStorage()->size());
    double execTime;
    std::cout << "The datavectors for training and the right hand side have been created" << std::endl;

    ARFFTool.readTrainingData(tfileTrain, training);
    std::cout << "The training vector has been initialized" << std::endl;
    ARFFTool.readClasses(tfileTrain, classes);
    std::cout << "The class training vector has been initialized" << std::endl;

    // init the Systemmatrix Functor
    OperationMatrix* C;
    if (this->StiffnessMode  == "L" )
    {
    	C = this->myGrid->createOperationLaplace();
    }
    else if (this->StiffnessMode == "I")
    {
    	C = this->myGrid->createOperationIdentity();
    }
    else
    {
    	C = this->myGrid->createOperationLaplace();
    }
    DMSystemMatrix DMMatrix(*this->myGrid, training, *C, this->lambda);
    std::cout << "Instance of the matrix functor has been created and initialized" << std::endl;
    // generate the rhs of the equation
    DMMatrix.generateb(classes, rhs);
    std::cout << "The rhs of the equation has been initialized" << std::endl;

    // get a CG
    ConjugateGradients myCG(this->IterationMax, this->epsilon);
    std::cout << "An instance of the CG method has been created" << std::endl;

    // slove the system of linear equations
    SGppStopwatch* myStopwatch = new SGppStopwatch();
    myStopwatch->start();
    myCG.solve(DMMatrix, alpha, rhs, false, false, -1.0);
    execTime = myStopwatch->stop();

    delete C;

    // Write the data of CG
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "Needed time: " << execTime << " seconds" << std::endl;
    std::cout << "Needed iterations: " << myCG.getNumberIterations() << std::endl;
    std::cout << "Final norm of residuum: " << myCG.getResiduum() << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << std::endl;
}

double Classifier::applyTestdata(DataVector& alpha, std::string tfileTest)
{
    double correct;

	ARFFTools ARFFTool;
	DataMatrix test(this->testinstancesNo, this->dim);
    DataVector testclasses(this->testinstancesNo);
    std::cout << "the test datavectors have been created" << std::endl;

    ARFFTool.readTrainingData(tfileTest, test);
    ARFFTool.readClasses(tfileTest, testclasses);
    std::cout << "the test datavectors have been initialized" << std::endl;

    std::cout << "start evaluating the test instances" << std::endl;

    OperationTest* myTest = myGrid->createOperationTest();
    correct = myTest->test(alpha, test, testclasses);
    delete myTest;

    std::cout << "finishes evaluating the test instances" << std::endl;

    std::cout << "Correctly classified elements: " << (correct/((double)test.getSize()))*100.0 << " %" << std::endl;

    return correct/((double)test.getSize());
}

}
