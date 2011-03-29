/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/pde/finance/OperationGammaLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
{

OperationGammaLinearBoundary::OperationGammaLinearBoundary(GridStorage* storage, DataMatrix& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLinearBoundary::~OperationGammaLinearBoundary()
{
}

void OperationGammaLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearBoundary func(this->storage);
	sweep<PhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearBoundary func(this->storage);
	sweep<PhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	XPhidPhiUpBBLinearBoundary func(this->storage);
	sweep<XPhidPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	XPhidPhiDownBBLinearBoundary func(this->storage);
	sweep<XPhidPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinearBoundary func(this->storage);
	sweep<XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinearBoundary func(this->storage);
	sweep<XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	SqXdPhidPhiUpBBLinearBoundary func(this->storage);
	sweep<SqXdPhidPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	SqXdPhidPhiDownBBLinearBoundary func(this->storage);
	sweep<SqXdPhidPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
