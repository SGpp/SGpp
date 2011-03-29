/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/LaserHeatEquationParabolicPDESolverSystemParallelOMP.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

LaserHeatEquationParabolicPDESolverSystemParallelOMP::LaserHeatEquationParabolicPDESolverSystemParallelOMP(Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode) : laser_velo_(1.0), HeatEquationParabolicPDESolverSystemParallelOMP(SparseGrid, alpha, a, TimestepSize, OperationMode)
{
}

LaserHeatEquationParabolicPDESolverSystemParallelOMP::~LaserHeatEquationParabolicPDESolverSystemParallelOMP()
{
}

void LaserHeatEquationParabolicPDESolverSystemParallelOMP::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void LaserHeatEquationParabolicPDESolverSystemParallelOMP::startTimestep()
{
}

}
