/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLDLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XPhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
{

OperationLDLinearBoundary::OperationLDLinearBoundary(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLDLinearBoundary::~OperationLDLinearBoundary()
{
}

void OperationLDLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// X * phi * phi
	XPhiPhiUpBBLinearBoundary func(this->storage);
	sweep<XPhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLDLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// X * phi * phi
	XPhiPhiDownBBLinearBoundary func(this->storage);
	sweep<XPhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
