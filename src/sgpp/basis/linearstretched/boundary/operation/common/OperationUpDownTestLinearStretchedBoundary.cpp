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
//	detail::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
//	detail::SqXdPhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::SqXdPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x * dphi * phi
	detail::XdPhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<detail::XdPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x * phi * dphi
//	detail::XPhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::XPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// dphi * phi
//	detail::DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// phi * dphi
//	detail::PhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::PhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationUpDownTestLinearStretchedBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
//	detail::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
//	detail::SqXdPhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::SqXdPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x * dphi * phi
	detail::XdPhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<detail::XdPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x * phi * dphi
//	detail::XPhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::XPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	//  dphi * phi
//	detail::DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	//  phi * dphi
//	detail::PhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sweep<detail::PhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
