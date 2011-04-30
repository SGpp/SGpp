/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretched.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationGammaLinearStretched::OperationGammaLinearStretched(GridStorage* storage, DataMatrix& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLinearStretched::~OperationGammaLinearStretched()
{
}

void OperationGammaLinearStretched::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretched func(this->storage);
	sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretched func(this->storage);
	sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	XPhidPhiUpBBLinearStretched func(this->storage);
	sweep<XPhidPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	XPhidPhiDownBBLinearStretched func(this->storage);
	sweep<XPhidPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinearStretched func(this->storage);
	sweep<XdPhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinearStretched func(this->storage);
	sweep<XdPhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	SqXdPhidPhiUpBBLinearStretched func(this->storage);
	sweep<SqXdPhidPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	SqXdPhidPhiDownBBLinearStretched func(this->storage);
	sweep<SqXdPhidPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
