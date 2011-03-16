/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLELinearBoundary.hpp"

//#include "basis/linear/boundary/algorithm_sweep/DPhidPhiDownBBLinearBoundary.hpp"
//#include "basis/linear/boundary/algorithm_sweep/DPhidPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/common/DowndPhidPhiBBIterativeLinearBoundary.hpp"
#include "basis/linear/boundary/common/UpdPhidPhiBBIterativeLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::base;

namespace sg
{

OperationLELinearBoundary::OperationLELinearBoundary(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLELinearBoundary::~OperationLELinearBoundary()
{
}

void OperationLELinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * dphi
	UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
	myUp(alpha, result, dim);
}

void OperationLELinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * dphi
	DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
	myDown(alpha, result, dim);
}

}
