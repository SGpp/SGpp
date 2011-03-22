/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "basis/linearstretched/boundary/operation/pde/OperationLaplaceLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/common/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/common/UpdPhidPhiBBIterativeLinearStretchedBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include "grid/common/Stretching.hpp"

namespace sg
{

OperationLaplaceLinearStretchedBoundary::OperationLaplaceLinearStretchedBoundary(GridStorage* storage) : UpDownOneOpDim(storage)
{
}

OperationLaplaceLinearStretchedBoundary::~OperationLaplaceLinearStretchedBoundary()
{
}

void OperationLaplaceLinearStretchedBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearStretchedBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearStretchedBoundary::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	DowndPhidPhiBBIterativeLinearStretchedBoundary myDown(this->storage);
	myDown(alpha, result, dim);
}

void OperationLaplaceLinearStretchedBoundary::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	UpdPhidPhiBBIterativeLinearStretchedBoundary myUp(this->storage);
	myUp(alpha, result, dim);
}

}
