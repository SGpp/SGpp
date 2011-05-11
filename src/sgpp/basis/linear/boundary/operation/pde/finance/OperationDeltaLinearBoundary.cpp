/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
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
	PhiPhiUpBBLinearBoundary func(this->storage);
	sweep<PhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearBoundary func(this->storage);
	sweep<PhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinearBoundary func(this->storage);
	sweep<XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinearBoundary func(this->storage);
	sweep<XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
