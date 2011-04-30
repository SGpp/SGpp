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
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationDeltaLinearStretched::OperationDeltaLinearStretched(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLinearStretched::~OperationDeltaLinearStretched()
{
}

void OperationDeltaLinearStretched::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretched func(this->storage);
	sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinearStretched::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretched func(this->storage);
	sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinearStretched::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinearStretched func(this->storage);
	sweep<XdPhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinearStretched::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinearStretched func(this->storage);
	sweep<XdPhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
