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

namespace sg
{
namespace pde
{

sg::pde::OperationUpDownTestLinearStretchedBoundary::OperationUpDownTestLinearStretchedBoundary(sg::base::GridStorage* storage) : storage(storage)
{
}

sg::pde::OperationUpDownTestLinearStretchedBoundary::~OperationUpDownTestLinearStretchedBoundary()
{
}

void sg::pde::OperationUpDownTestLinearStretchedBoundary::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	this->updown(alpha, result);
}

void sg::pde::OperationUpDownTestLinearStretchedBoundary::updown(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector beta(result.getSize());

	this->updown(alpha, beta, storage->dim() - 1);

	result.add(beta);
}

void sg::pde::OperationUpDownTestLinearStretchedBoundary::updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		sg::base::DataVector temp(alpha.getSize());
		up(alpha, temp, dim);
		updown(temp, result, dim-1);


		// Same from the other direction
		sg::base::DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1);
		down(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		up(alpha, result, dim);

		sg::base::DataVector temp(alpha.getSize());
		down(alpha, temp, dim);

		result.add(temp);
	}
}

void sg::pde::OperationUpDownTestLinearStretchedBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
//	sg::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
//	sg::finance::SqXdPhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::SqXdPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x * dphi * phi
	sg::finance::XdPhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sg::base::sweep<sg::finance::XdPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// x * phi * dphi
//	sg::finance::XPhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::XPhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// dphi * phi
//	sg::finance::DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	// phi * dphi
//	sg::finance::PhidPhiUpBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::PhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void sg::pde::OperationUpDownTestLinearStretchedBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
//	sg::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
//	sg::finance::SqXdPhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::SqXdPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x * dphi * phi
	sg::finance::XdPhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sg::base::sweep<sg::finance::XdPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	// x * phi * dphi
//	sg::finance::XPhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::XPhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	//  dphi * phi
//	sg::finance::DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	//  phi * dphi
//	sg::finance::PhidPhiDownBBLinearStretchedBoundary func(this->storage);
//	sg::base::sweep<sg::finance::PhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
