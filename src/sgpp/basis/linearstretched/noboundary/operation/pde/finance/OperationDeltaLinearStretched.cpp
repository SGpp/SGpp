/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiUpBBLinearStretched.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{
namespace finance
{

OperationDeltaLinearStretched::OperationDeltaLinearStretched(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLinearStretched::~OperationDeltaLinearStretched()
{
}

void OperationDeltaLinearStretched::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	sg::pde::PhiPhiUpBBLinearStretched func(this->storage);
	sg::base::sweep<sg::pde::PhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinearStretched::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	sg::pde::PhiPhiDownBBLinearStretched func(this->storage);
	sg::base::sweep<sg::pde::PhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinearStretched::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinearStretched func(this->storage);
	sg::base::sweep<XdPhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinearStretched::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinearStretched func(this->storage);
	sg::base::sweep<XdPhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
