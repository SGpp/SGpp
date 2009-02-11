/*****************************************************************************/
/* This file is part of sg++, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sg++ is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sg++ is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sg++; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef OPERATIONLAPLACEMODLINEAR_HPP
#define OPERATIONLAPLACEMODLINEAR_HPP

#include "basis/modlinear/algorithm_sweep/LaplaceDownGradientModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/LaplaceDownModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/LaplaceUpGradientModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/LaplaceUpModLinear.hpp"

#include "operation/OperationMatrix.hpp"

#include "algorithm/UnidirGradient.hpp"
#include "algorithm/sweep.hpp"

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

class OperationLaplaceModLinear : public UnidirGradient, public OperationMatrix
{
public:
	OperationLaplaceModLinear(GridStorage* storage) : UnidirGradient(storage) {}
	virtual ~OperationLaplaceModLinear() {}

	virtual void mult(DataVector& alpha, DataVector& result)
	{
		this->updown(alpha, result);
	}

protected:

	virtual void up(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceUpModLinear func(this->storage);
		sweep<detail::LaplaceUpModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

	virtual void down(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceDownModLinear func(this->storage);
		sweep<detail::LaplaceDownModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

	virtual void downGradient(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceDownGradientModLinear func(this->storage);
		sweep<detail::LaplaceDownGradientModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

	virtual void upGradient(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceUpGradientModLinear func(this->storage);
		sweep<detail::LaplaceUpGradientModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}
};

}

#endif /* OPERATIONLAPLACEMODLINEAR_HPP */
