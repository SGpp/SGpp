/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{
namespace finance
{

OperationDeltaLinear::OperationDeltaLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLinear::~OperationDeltaLinear()
{
}

void OperationDeltaLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	sg::pde::PhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	sg::pde::PhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
