/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/LaserHeatEquationParabolicPDESolverSystemParallelOMP.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "tools/common/StdNormalDistribution.hpp"

#include <string>
#include <sstream>
#include <cmath>

#ifndef PI
#define PI 3.14159265
#endif

namespace sg
{

LaserHeatEquationParabolicPDESolverSystemParallelOMP::LaserHeatEquationParabolicPDESolverSystemParallelOMP(double beam_velocity, double heat_sigma, size_t max_level, Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode) : beam_velocity_(beam_velocity), heat_sigma_(heat_sigma), max_level_(max_level), laser_x_offset_(0.0), laser_x_start_(0.25), laser_x_last_(0.25), HeatEquationParabolicPDESolverSystemParallelOMP(SparseGrid, alpha, a, TimestepSize, OperationMode)
{
}

LaserHeatEquationParabolicPDESolverSystemParallelOMP::~LaserHeatEquationParabolicPDESolverSystemParallelOMP()
{
}

void LaserHeatEquationParabolicPDESolverSystemParallelOMP::finishTimestep(bool isLastTimestep)
{
	double heat = 4.0;
	double xadd = 0.0;
	bool y_add = true;
	bool x_add = true;

	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

	// apply new laser position
	StdNormalDistribution myNormDistr;
	OperationHierarchisation* myHierarchisation = this->BoundGrid->createOperationHierarchisation();
	myHierarchisation->doDehierarchisation(*this->alpha_complete);

	DataVector laser_update(this->BoundGrid->getStorage()->size());
	BoundingBox* myBoundingBox = new BoundingBox(*(this->BoundGrid->getStorage()->getBoundingBox()));

	this->laser_x_offset_ += (this->beam_velocity_*this->TimestepSize);

	int int_x_offset = (int)laser_x_offset_;
	if ((this->laser_x_offset_ - ((double)int_x_offset)) <= 0.25)
	{
		y_add = true;
		x_add = true;
		xadd = 0.0;
	}
	if ((this->laser_x_offset_ - ((double)int_x_offset)) > 0.25)
	{
		y_add = true;
		x_add = false;
		xadd = 0.5;
	}
	if ((this->laser_x_offset_ - ((double)int_x_offset)) > 0.5)
	{
		y_add = false;
		x_add = true;
		xadd = 0.5;
	}
	if ((this->laser_x_offset_ - ((double)int_x_offset)) > 0.75)
	{
		y_add = false;
		x_add = false;
		xadd = 0.0;
	}

	double pos_x = 0.0;
	if (x_add == true)
	{
		pos_x = this->laser_x_start_ + xadd + (0.25*sin(this->laser_x_offset_*2.0*PI));
	}
	else
	{
		pos_x = this->laser_x_start_ + xadd - (0.25*sin(this->laser_x_offset_*2.0*PI));
	}
	double radi = (0.25*0.25)-((pos_x-0.5)*(pos_x-0.5));
	if (radi < 0.0)
	{
		radi = 0.0;
	}
	double pos_y = sqrt(radi);

	if (y_add == true)
	{
		pos_y = pos_y + 0.5;
	}
	else
	{
		pos_y = 0.5 - pos_y;
	}

	//std::cout << std::endl << std::endl << pos_x << " " << pos_y << std::endl << std::endl;

	double* dblFuncValues = new double[2];
	for (size_t i = 0; i < this->BoundGrid->getStorage()->size(); i++)
	{
		std::string coords = this->BoundGrid->getStorage()->get(i)->getCoordsStringBB(*(myBoundingBox));
		std::stringstream coordsStream(coords);

		for (size_t j = 0; j < 2; j++)
		{
			coordsStream >> dblFuncValues[j];
		}

		// check if coordinates at starting point of laser
		laser_update[i] =  heat*(myNormDistr.getDensity(dblFuncValues[0], pos_x, this->heat_sigma_)*myNormDistr.getDensity(dblFuncValues[1], pos_y, this->heat_sigma_));

		//boundaries are set to zero
		if (dblFuncValues[0] == 0.0 || dblFuncValues[1] == 0.0)
		{
			laser_update[i] = 0.0;
		}
	}
	delete[] dblFuncValues;
	delete myBoundingBox;

	this->laser_x_last_ = pos_x;

	// combine last solution and laser update
	for (size_t i = 0; i < this->BoundGrid->getStorage()->size(); i++)
	{
		this->alpha_complete->set(i, std::max< double >(this->alpha_complete->get(i), laser_update[i]));
	}

	// do hierarchisation
	myHierarchisation->doHierarchisation(*this->alpha_complete);
	delete myHierarchisation;

	///////////////////////////////////////////////////
	// Start integrated refinement & coarsening
	///////////////////////////////////////////////////

	size_t originalGridSize = this->BoundGrid->getStorage()->size();

	// Coarsen the grid
	GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

	size_t numRefines = myGenerator->getNumberOfRefinablePoints();
	SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(this->alpha_complete, numRefines, 0.01);
	myGenerator->refineMaxLevel(myRefineFunc, this->max_level_);
	this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
	delete myRefineFunc;

	if (this->BoundGrid->getStorage()->getNumInnerPoints() > 100)
	{
		size_t numCoarsen = myGenerator->getNumberOfRemoveablePoints();
		SurplusCoarseningFunctor* myCoarsenFunctor = new SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, 0.001);
		myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
		delete myCoarsenFunctor;
	}

	delete myGenerator;

	///////////////////////////////////////////////////
	// End integrated refinement & coarsening
	///////////////////////////////////////////////////

	// rebuild the inner grid + coefficients
	this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
}

void LaserHeatEquationParabolicPDESolverSystemParallelOMP::startTimestep()
{
}

}
