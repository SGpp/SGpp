/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationDeltaLinearStretchedBoundary::OperationDeltaLinearStretchedBoundary(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLinearStretchedBoundary::~OperationDeltaLinearStretchedBoundary()
{
}

void OperationDeltaLinearStretchedBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearStretchedBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearStretchedBoundary::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<XdPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearStretchedBoundary::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<XdPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
