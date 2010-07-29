/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
	size_t levels = 3;
	size_t dim = 1;
	size_t numGridPoints;

	double coarsenThreshold = 0.01;
	size_t numCoarsen = 2;

	sg::Grid* myGrid;

	std::cout << std::endl;
	std::cout << "Starting the Refine/Coarsen Test" << std::endl;
	std::cout << "================================" << std::endl;
	std::cout << std::endl;
	std::cout << "levels:" << levels << std::endl;
	std::cout << "dim:   " << dim << std::endl;
	std::cout << std::endl;

	std::cout << "start constructing regular grid" << std::endl;
	myGrid = new sg::LinearTrapezoidBoundaryGrid(dim);
	std::cout << "A LinearTrapezoidBoundaryGrid was created" << std::endl;

	sg::GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(levels);
	delete myGenerator;
	std::cout << levels << " levels were added to the above created grid" << std::endl;

	numGridPoints = myGrid->getStorage()->size();
	std::cout << "the grid has " << numGridPoints << " gridpoints" << std::endl;
	std::cout << "finished construction regular grid" << std::endl << std::endl;

	std::cout << "===========================================================" << std::endl;
	std::cout << "Coarsening parameters are as followed:" << std::endl;
	std::cout << "threshold: " << coarsenThreshold << std::endl;
	std::cout << "No. Coarsenings: " << numCoarsen <<  std::endl;
	std::cout << "===========================================================" << std::endl;
	std::cout << std::endl;

	std::string ser;
	myGrid->serialize(ser);
	std::cout << "the grid is:" << std::endl;
	std::cout << ser << std::endl;

	DataVector alpha(numGridPoints);
	alpha.setAll(1.0);
	alpha.set(5, 0.00001);
	alpha.set(6, 0.1);

	std::cout << "Orig alpha:" << std::endl;
	std::cout << alpha.toString() << std::endl << std::endl;

	sg::GridGenerator* myGeneratorCoarsen = myGrid->createGridGenerator();
	sg::SurplusCoarseningFunctor* myCoarsenFunctor = new sg::SurplusCoarseningFunctor(&alpha, numCoarsen, coarsenThreshold);

	myGeneratorCoarsen->coarsen(myCoarsenFunctor, &alpha);

	delete myGeneratorCoarsen;
	delete myCoarsenFunctor;

	ser = "";
	myGrid->serialize(ser);
	std::cout << "the new grid is:" << std::endl;
	std::cout << ser << std::endl;

	std::cout << "New alpha:" << std::endl;
	std::cout << alpha.toString() << std::endl << std::endl;

	return 0;
}
