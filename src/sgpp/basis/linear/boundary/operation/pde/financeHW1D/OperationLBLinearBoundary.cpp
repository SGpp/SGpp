/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
{

OperationLBLinearBoundary::OperationLBLinearBoundary(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLBLinearBoundary::~OperationLBLinearBoundary()
{
}

void OperationLBLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	DPhiPhiUpBBLinearBoundary func(this->storage);
	sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLBLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	DPhiPhiDownBBLinearBoundary func(this->storage);
	sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
