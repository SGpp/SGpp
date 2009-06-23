/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "application/classification/Classifier.hpp"
#include "exception/operation_exception.hpp"
#include "solver/cg/ApplyMatrix/classification/ApplyDMMatrix.hpp"
#include "solver/cg/ConjugateGradients.hpp"
#include "tools/classification/ARFFTools.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"
#include "grid/type/LinearBoundaryUScaledGrid.hpp"
#include "grid/type/ModLinearGrid.hpp"
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
		myGrid = new LinearBoundaryUScaledGrid(this->dim);
		std::cout << "A LinearBoundaryUScaledGrid was created" << std::endl;
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

	myGrid->createGridGenerator()->regular(this->levels);
	std::cout << levels << " levels were added to the above created grid" << std::endl;
}

void Classifier::trainGrid(DataVector& alpha, std::string tfileTrain)
{
    ARFFTools ARFFTool;
	DataVector training(this->instancesNo, this->dim);
    DataVector classes(this->instancesNo);
    DataVector rhs(this->myGrid->getStorage()->size());
    std::cout << "The datavectors for training and the right hand side have been created" << std::endl;

    ARFFTool.readTrainingData(tfileTrain, training);
    std::cout << "The training vector has been initialized" << std::endl;
    ARFFTool.readClasses(tfileTrain, classes);
    std::cout << "The class training vector has been initialized" << std::endl;

    // init the Systemmatrix Functor
    ApplyDMMatrix DMMatrix(this->myGrid, this->StiffnessMode, this->lambda);
    std::cout << "Instance of the matrix functor has been created and initialized" << std::endl;
    // generate the rhs of the equation
    DMMatrix.generateb(training, classes, rhs);
    std::cout << "The rhs of the equation has been initialized" << std::endl;

    // get a CG
    ConjugateGradients<ApplyDMMatrix> myCG(this->IterationMax, this->epsilon);
    std::cout << "An instance of the CG method has been created" << std::endl;

    // slove the system of linear equations
    myCG.solve(DMMatrix, alpha, training, rhs, true, false);

    // Write the data of CG
    std::cout << "Needed iterations: " << myCG.getNumberIterations() << std::endl;
}

double Classifier::applyTestdata(DataVector& alpha, std::string tfileTest)
{
    double correct;

	ARFFTools ARFFTool;
	DataVector test(this->testinstancesNo, this->dim);
    DataVector testclasses(this->testinstancesNo);
    std::cout << "the test datavectors have been created" << std::endl;

    ARFFTool.readTrainingData(tfileTest, test);
    ARFFTool.readClasses(tfileTest, testclasses);
    std::cout << "the test datavectors have been initialized" << std::endl;

    std::cout << "start evaluating the test instances" << std::endl;
    correct = myGrid->createOperationEval()->test(alpha, test, testclasses);
    std::cout << "finishes evaluating the test instances" << std::endl;

    std::cout << "Correctly classified elements: " << (correct/((double)test.getSize()))*100.0 << " percentage" << std::endl;

    return correct/((double)test.getSize());
}

}
