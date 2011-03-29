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

#ifndef PI
#define PI 3.14159265
#endif

namespace sg
{

LaserHeatEquationParabolicPDESolverSystemParallelOMP::LaserHeatEquationParabolicPDESolverSystemParallelOMP(double beam_velocity, double heat_sigma, size_t max_level, Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode) : beam_velocity_(beam_velocity), heat_sigma_(heat_sigma), max_level_(max_level), laser_x_offset_(0.0), laser_x_start_(0.25), laser_x_last_(0.25), y_sign_switch(0), HeatEquationParabolicPDESolverSystemParallelOMP(SparseGrid, alpha, a, TimestepSize, OperationMode)
{
}

LaserHeatEquationParabolicPDESolverSystemParallelOMP::~LaserHeatEquationParabolicPDESolverSystemParallelOMP()
{
}

void LaserHeatEquationParabolicPDESolverSystemParallelOMP::finishTimestep(bool isLastTimestep)
{
	double heat = 4.0;

	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

	// apply new laser position
	StdNormalDistribution myNormDistr;
	OperationHierarchisation* myHierarchisation = this->BoundGrid->createOperationHierarchisation();
	myHierarchisation->doDehierarchisation(*this->alpha_complete);

//	DataVector laser_update(this->alpha_complete->getSize());
//
//	this->laser_x_offset_ += (this->beam_velocity_*this->TimestepSize);
//	double pos_x = this->laser_x_start_ + sin(this->laser_x_offset_*2.0*PI);
//	if ((pos_x > 0.5) && (this->laser_x_last_ <= 0.5))
//	{
//		this->y_sign_switch++;
//	}
//	if ((pos_x < 0.5) && (this->laser_x_last_ >= 0.5))
//	{
//		this->y_sign_switch++;
//	}
//	double pos_y = sqrt((0.25*0.25)-((pos_x-0,5)*(pos_x-0,5)));
//	pos_y = pos_y*(pow(-1.0, (double)this->y_sign_switch));
//	pos_y += 0.5;
//
//	double* dblFuncValues = new double[2];
//	for (size_t i = 0; i < this->BoundGrid->getStorage()->size(); i++)
//	{
//		std::string coords = this->BoundGrid->getStorage()->get(i)->getCoordsStringBB(*(this->BoundGrid->getStorage()->getBoundingBox()));
//		std::stringstream coordsStream(coords);
//
//		for (size_t j = 0; j < 2; j++)
//		{
//			coordsStream >> dblFuncValues[j];
//		}
//
//		// check if coordinates at starting point of laser
//		laser_update[i] =  heat*(myNormDistr.getDensity(dblFuncValues[0], pos_x, this->heat_sigma_)*myNormDistr.getDensity(dblFuncValues[1], pos_y, this->heat_sigma_));
//
//		//boundaries are set to zero
//		if (dblFuncValues[0] == 0.0 || dblFuncValues[1] == 0.0)
//		{
//			laser_update[i] = 0.0;
//		}
//	}
//	delete[] dblFuncValues;
//
//	this->laser_x_last_ = pos_x;
//
//	// combine last solution and laser update
//	for (size_t i = 0; this->alpha_complete->getSize(); i++)
//	{
//		this->alpha_complete->set(i, std::max< double >(this->alpha_complete->get(i), laser_update[i]));
//	}

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
	SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(this->alpha_complete, numRefines, 0.0001);
	myGenerator->refineMaxLevel(myRefineFunc, this->max_level_);
	this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
	delete myRefineFunc;

	if (this->BoundGrid->getStorage()->getNumInnerPoints() > 100)
	{
		size_t numCoarsen = myGenerator->getNumberOfRemoveablePoints();
		SurplusCoarseningFunctor* myCoarsenFunctor = new SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, 0.00001);
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
