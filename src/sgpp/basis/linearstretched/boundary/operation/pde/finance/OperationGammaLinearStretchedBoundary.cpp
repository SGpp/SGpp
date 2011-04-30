/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/XPhidPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/XPhidPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretchedBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationGammaLinearStretchedBoundary::OperationGammaLinearStretchedBoundary(GridStorage* storage, DataMatrix& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLinearStretchedBoundary::~OperationGammaLinearStretchedBoundary()
{
}

void OperationGammaLinearStretchedBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	XPhidPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<XPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	XPhidPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<XPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<XdPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<XdPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	SqXdPhidPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<SqXdPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	SqXdPhidPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<SqXdPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
