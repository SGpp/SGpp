/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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

#include "basis/linear/modlinear/operation/pde/OperationLaplaceModLinear.hpp"

#include "basis/linear/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp"
#include "basis/linear/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp"
#include "basis/linear/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp"
#include "basis/linear/modlinear/algorithm_sweep/PhiPhiUpModLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLaplaceModLinear::OperationLaplaceModLinear(GridStorage* storage) : UnidirGradient(storage)
{
}

OperationLaplaceModLinear::~OperationLaplaceModLinear()
{
}

void OperationLaplaceModLinear::mult(DataVector& alpha, DataVector& result)
{
	this->updown(alpha, result);
}

void OperationLaplaceModLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	detail::PhiPhiUpModLinear func(this->storage);
	sweep<detail::PhiPhiUpModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	detail::PhiPhiDownModLinear func(this->storage);
	sweep<detail::PhiPhiDownModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	detail::dPhidPhiDownModLinear func(this->storage);
	sweep<detail::dPhidPhiDownModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	detail::dPhidPhiUpModLinear func(this->storage);
	sweep<detail::dPhidPhiUpModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

}

