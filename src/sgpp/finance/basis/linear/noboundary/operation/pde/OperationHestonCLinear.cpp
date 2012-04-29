/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonCLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg
{
namespace finance
{

OperationHestonCLinear::OperationHestonCLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationHestonCLinear::~OperationHestonCLinear()
{
}

void OperationHestonCLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonCLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonCLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonCLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
