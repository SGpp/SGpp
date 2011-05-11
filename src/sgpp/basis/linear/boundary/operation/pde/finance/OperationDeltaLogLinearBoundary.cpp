/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLogLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
{

OperationDeltaLogLinearBoundary::OperationDeltaLogLinearBoundary(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLogLinearBoundary::~OperationDeltaLogLinearBoundary()
{
}

void OperationDeltaLogLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearBoundary func(this->storage);
	sweep<PhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearBoundary func(this->storage);
	sweep<PhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearBoundary::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinearBoundary func(this->storage);
	sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearBoundary::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinearBoundary func(this->storage);
	sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
