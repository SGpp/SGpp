/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "basis/linearstretched/boundary/operation/common/OperationUpDownTestLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/XPhidPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/XPhidPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/DPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/DPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhidPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhidPhiUpBBLinearStretchedBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::finance;

namespace sg
{
namespace pde
{

OperationUpDownTestLinearStretchedBoundary::OperationUpDownTestLinearStretchedBoundary(GridStorage* storage) : storage(storage)
{
}

OperationUpDownTestLinearStretchedBoundary::~OperationUpDownTestLinearStretchedBoundary()
{
}

void OperationUpDownTestLinearStretchedBoundary::mult(DataVector& alpha, DataVector& result)
{
	this->updown(alpha, result);
}

void OperationUpDownTestLinearStretchedBoundary::updown(DataVector& alpha, DataVector& result)
{
	DataVector beta(result.getSize());

	this->updown(alpha, beta, storage->dim() - 1);

	result.add(beta);
}

void OperationUpDownTestLinearStretchedBoundary::updown(DataVector& alpha, DataVector& result, size_t dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		up(alpha, temp, dim);
		updown(temp, result, dim-1);


		// Same from the other direction
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1);
		down(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		up(alpha, result, dim);

		DataVector temp(alpha.getSize());
		down(alpha, temp, dim);

		result.add(temp);
	}
}

void OperationUpDownTestLinearStretchedBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
//	PhiPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
//	SqXdPhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<SqXdPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x * dphi * phi
	XdPhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<XdPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x * phi * dphi
//	XPhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<XPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// dphi * phi
//	DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// phi * dphi
//	PhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<PhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationUpDownTestLinearStretchedBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
//	PhiPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
//	SqXdPhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<SqXdPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x * dphi * phi
	XdPhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<XdPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x * phi * dphi
//	XPhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<XPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	//  dphi * phi
//	DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	//  phi * dphi
//	PhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<PhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
