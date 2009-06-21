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
	this->lambda = 0.0;
	this->StiffnessMode = "L";
	this->IterationMax = 1000;
	this->epsilon = 0.0;
}

Classifier::~Classifier()
{
	delete myGrid;
}

void Classifier::trainNtestRegular(std::string tfileTrain, std::string tfileTest, size_t level, double lambda, char GridType, char StiffMode, double epsilon, size_t imax)
{
	ARFFTools ARFFTool;

	this->levels = level;
	this->GridType = GridType;
	this->lambda = lambda;
	this->StiffnessMode = StiffMode;
	this->IterationMax = imax;
	this->epsilon = epsilon;
	this->dim = ARFFTool.getDimension(tfileTrain);
	this->instancesNo = ARFFTool.getDimension(tfileTrain);

	// create the grid
	createRegularGrid();

	DataVector alpha(this->myGrid->getStorage()->size());
	alpha.set(0.0);

	// start the train process
	trainGrid(alpha, tfileTrain);

	// start test process
	applyTestdata(alpha, tfileTest);
}

void Classifier::createRegularGrid()
{
	switch (this->GridType)
	{
	case "N":
		myGrid = new LinearGrid(this->dim);
		break;
	case "U":
		myGrid = new LinearBoundaryUScaledGrid(this->dim);
		break;
	case "B":
		myGrid = new LinearBoundaryGrid(this->dim);
		break;
	case "E":
		myGrid = new ModLinearGrid(this->dim);
		break;
	default:
		throw new operation_exception("No valid GridType was specified!");
		break;
	}

	myGrid->createGridGenerator()->regular(this->levels);
}

void Classifier::trainGrid(DataVector& alpha, std::string tfileTrain)
{
    ARFFTools ARFFTool;
	DataVector training(this->instancesNo, this->dim);
    DataVector classes(this->instancesNo);
    DataVector rhs(this->myGrid->getStorage()->size());

    ARFFTool.readTrainingData(tfileTrain, training);
    ARFFTool.readClasses(tfileTrain, classes);

    // init the Systemmatrix Functor
    ApplyDMMatrix DMMatrix(this->myGrid, this->StiffnessMode, this->lambda);
    // generate the rhs of the equation
    DMMatrix.generateb(training, classes, rhs);

    // get a CG
    ConjugateGradients<ApplyDMMatrix> myCG(this->IterationsMax, this->epsilon);

    // slove the system of linear equations
    myCG.solve(DDMatrix, alpha, training, rhs, false, true);

    // Write the data of CG
    std::cout << "Final norm of residual: " << myCG.getFinalResiduum() << std::endl;
    std::cout << "Needed iterations: " << myCG.getNumberIterations() << std::endl;
}

double Classifier::applyTestdata(DataVector& alpha, std::string tfileTest)
{
    double correct;

	ARFFTools ARFFTool;
	DataVector test(this->instancesNo, this->dim);
    DataVector testclasses(this->instancesNo);

    ARFFTool.readTrainingData(tfileTest, test);
    ARFFTool.readClasses(tfileTest, testclasses);

    correct = myGrid->createOperationEval()->test(alpha, test, testclasses);

    std::cout << "Correctly classified elements: " << correct/((double)test.getSize()) << " percentage" << std::endl;

    return correct/((double)test.getSize());
}

}
