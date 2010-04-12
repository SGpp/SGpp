/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>

namespace sg
{

OperationGammaLinear::OperationGammaLinear(GridStorage* storage, DataVector& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLinear::~OperationGammaLinear()
{
}

void OperationGammaLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinear func(this->storage);
	sweep<detail::PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinear func(this->storage);
	sweep<detail::PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiUpBBLinear func(this->storage);
	sweep<detail::XPhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiDownBBLinear func(this->storage);
	sweep<detail::XPhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinear::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiUpBBLinear func(this->storage);
	sweep<detail::XdPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinear::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiDownBBLinear func(this->storage);
	sweep<detail::XdPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiUpBBLinear func(this->storage);
	sweep<detail::SqXdPhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiDownBBLinear func(this->storage);
	sweep<detail::SqXdPhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
