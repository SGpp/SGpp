/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationDeltaLinearBoundary::OperationDeltaLinearBoundary(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLinearBoundary::~OperationDeltaLinearBoundary()
{
}

void OperationDeltaLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::XPhidPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::XPhidPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
