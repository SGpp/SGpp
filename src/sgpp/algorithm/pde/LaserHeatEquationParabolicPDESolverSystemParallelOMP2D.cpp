/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/LaserHeatEquationParabolicPDESolverSystemParallelOMP2D.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "tools/common/StdNormalDistribution.hpp"
#include "base/operation/BaseOpFactory.hpp"

#include <string>
#include <sstream>
#include <cmath>

#ifndef PI
#define PI 3.14159265
#endif

namespace sg
{
namespace pde
{

LaserHeatEquationParabolicPDESolverSystemParallelOMP2D::LaserHeatEquationParabolicPDESolverSystemParallelOMP2D(double beam_velocity, double heat_sigma, size_t max_level, double heat, double refine_threshold, double coarsen_threshold, sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, double a, double TimestepSize, std::string OperationMode) : beam_velocity_(beam_velocity), heat_sigma_(heat_sigma), max_level_(max_level), heat_(heat), refine_threshold_(refine_threshold), coarsen_threshold_(coarsen_threshold), done_steps_(0), HeatEquationParabolicPDESolverSystemParallelOMP(SparseGrid, alpha, a, TimestepSize, OperationMode)
{
}

LaserHeatEquationParabolicPDESolverSystemParallelOMP2D::~LaserHeatEquationParabolicPDESolverSystemParallelOMP2D()
{
}

void LaserHeatEquationParabolicPDESolverSystemParallelOMP2D::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

	// apply new laser position
	base::StdNormalDistribution myNormDistr;
	base::OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->BoundGrid);
	myHierarchisation->doDehierarchisation(*this->alpha_complete);

	base::DataVector laser_update(this->BoundGrid->getStorage()->size());
	base::BoundingBox* myBoundingBox = new base::BoundingBox(*(this->BoundGrid->getStorage()->getBoundingBox()));

	this->done_steps_++;
	double angle = ((this->beam_velocity_*(double(this->done_steps_)*this->TimestepSize)) * 360) - 180;

	double pos_x = ((cos((angle * PI) / 180.0))*0.25)+0.5;
	double pos_y = ((sin((angle * PI) / 180.0))*0.25)+0.5;

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
		laser_update[i] =  this->heat_*(myNormDistr.getDensity(dblFuncValues[0], pos_x, this->heat_sigma_)*myNormDistr.getDensity(dblFuncValues[1], pos_y, this->heat_sigma_));

		//boundaries are set to zero
		if (dblFuncValues[0] == 0.0 || dblFuncValues[1] == 0.0)
		{
			laser_update[i] = 0.0;
		}
	}
	delete[] dblFuncValues;
	delete myBoundingBox;

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
	base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

	size_t numRefines = myGenerator->getNumberOfRefinablePoints();
	base::SurplusRefinementFunctor* myRefineFunc = new base::SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refine_threshold_);
	myGenerator->refineMaxLevel(myRefineFunc, this->max_level_);
	this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
	delete myRefineFunc;

	if (this->BoundGrid->getStorage()->getNumInnerPoints() > 100)
	{
		size_t numCoarsen = myGenerator->getNumberOfRemoveablePoints();
		base::SurplusCoarseningFunctor* myCoarsenFunctor = new base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsen_threshold_);
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

void LaserHeatEquationParabolicPDESolverSystemParallelOMP2D::startTimestep()
{
}

}
}
