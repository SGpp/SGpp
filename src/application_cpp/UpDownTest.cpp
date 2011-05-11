/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "exception/operation_exception.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"
#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/type/ModLinearGrid.hpp"
#include "basis/operations_factory.hpp"

#include <iostream>
#include <string>
#include <fstream>

#define PRINTMATRIX
#define PRINTDIFF

int main(int argc, char *argv[])
{
	size_t levels = 4;
	size_t dim = 1;
	size_t numGridPoints;

	sg::base::Grid* myGrid;
	sg::base::OperationMatrix* myUpDown;

	std::cout << std::endl;
	std::cout << "Starting the Up / Down Test" << std::endl;
	std::cout << "===========================" << std::endl;
	std::cout << std::endl;
	std::cout << "levels:" << levels << std::endl;
	std::cout << "dim:   " << dim << std::endl;
	std::cout << std::endl;

	std::cout << "start constructing regular grid" << std::endl;
	myGrid = new sg::base::LinearTrapezoidBoundaryGrid(dim);
	std::cout << "A LinearTrapezoidBoundaryGrid was created" << std::endl;

	sg::base::GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(levels);
	delete myGenerator;
	std::cout << levels << " levels were added to the above created grid" << std::endl;

	numGridPoints = myGrid->getStorage()->size();
	std::cout << "the grid has " << numGridPoints << " gridpoints" << std::endl;
	std::cout << "finished construction regular grid" << std::endl;
	//std::cout << "the grid is:" << std::endl;

	std::string ser;
	myGrid->serialize(ser);
	//std::cout << ser << std::endl;

	sg::base::BoundingBox* myBoundingBox;
	myBoundingBox = myGrid->getBoundingBox();
	sg::base::DimensionBoundary myInterval;

	myInterval.leftBoundary = 0.0;
	myInterval.rightBoundary = 100.0;
	myInterval.bDirichletLeft = false;
	myInterval.bDirichletRight = false;

	myBoundingBox->setBoundary(0, myInterval);

	std::cout << "Changed Bounding Box to: " << myInterval.leftBoundary << " to " << myInterval.rightBoundary << std::endl;

	DataVector alpha(numGridPoints);
	alpha.setAll(0.0);

	DataVector result(numGridPoints);
	result.setAll(0.0);

	DataMatrix UpDownMatrix(numGridPoints, numGridPoints);
	UpDownMatrix.setAll(0.0);

	myUpDown = sg::GridOperationFactory::createOperationUpDownTest(*myGrid);

	std::cout << "start constructing the operator's matrix" << std::endl;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		// init alpha
		alpha.setAll(0.0);
		alpha.set(i, 1.0);
		// init result
		result.setAll(0.0);

		myUpDown->mult(alpha, result);

		// copy data to opartor's matrix
		UpDownMatrix.setColumn(i, result);
	}
	std::cout << "finished constructing the operator's matrix" << std::endl;

#ifdef PRINTMATRIX
	std::cout << "The operator's matrix is:" << std::endl << std::endl;
	std::ofstream outfile;

	outfile.open("updown.test");

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			std::cout << std::scientific << UpDownMatrix.get(i, j) << " ";
			outfile << std::scientific << UpDownMatrix.get(i, j) << " ";
		}
		std::cout << std::endl;
		outfile << std::endl;
	}

	outfile.close();
#endif

#ifdef PRINTDIFF
	std::string file = "XdPhiPhi1D.txt";
	std::ifstream infile;

	infile.open(file.c_str());

	double filedata = 0.0;

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			infile >> filedata;
			std::cout << (UpDownMatrix.get(i, j) - filedata) << " ";
		}
		std::cout << std::endl;
	}

	infile.close();
#endif

	std::cout << std::endl;

	std::cout << "Test symmetry:" << std::endl;
	double symTest = 0.0;
	double tempSym = 0.0;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			tempSym = fabs(UpDownMatrix.get(i,j)-UpDownMatrix.get(j,i));
			if (tempSym > symTest)
			{
				symTest = tempSym;
			}
		}
	}
	std::cout << "Maximum symmetry error: " << symTest << std::endl << std::endl;

	std::cout << "Do a Test multiplication:" << std::endl << std::endl;

	alpha.setAll(0.0);
	alpha.set(1, 35.0);
	alpha.set(2, -17.5);
	alpha.set(4, -7.5);
	alpha.set(7, -5);
	alpha.set(14, -1.25);

	result.setAll(0.0);

	myUpDown->mult(alpha, result);

	std::cout << result.toString() << std::endl;

	std::cout << std::endl;

	/*sg::IOToolBonnSG* myImporter = new sg::IOToolBonnSG();

	DataVector serAlpha(0);
	std::string serGrid;
	bool hier;

	myImporter->readFile("sparse_bonn.grid", serGrid, serAlpha, hier);

	myImporter->writeFile("sgpp_bonn.grid", *myGrid, alpha, true);

	delete myImporter;*/
	delete myUpDown;
	delete myGrid;

#ifdef WINDOWS
	system("pause");
#endif
	return 0;
}
