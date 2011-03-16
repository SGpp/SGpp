/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef OPERATIONLAPLACEPREWAVELET_HPP
#define OPERATIONLAPLACEPREWAVELET_HPP

#include "basis/prewavelet/algorithm_sweep/LaplaceDownGradientPrewavelet.hpp"
#include "basis/prewavelet/algorithm_sweep/LaplaceUpGradientPrewavelet.hpp"
#include "basis/prewavelet/algorithm_sweep/LaplaceUpPrewavelet.hpp"
#include "algorithm/pde/UpDownOneOpDimWithShadow.hpp"

#include "operation/common/OperationMatrix.hpp"

#include "algorithm/common/sweep.hpp"

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"
using namespace sg::base;

#include<iostream>

namespace sg
{

/**
 * Implementation for linear functions of Laplace Operation, prewavelet grids without boundaries.
 * With prewavelets the calculation of the gradient part of the up down algorithm is the more complicated
 * one whereas the normal part is eased. For details on the implementation please refer to the documentation
 * of the detail-classes LaplaceDownGradientPrewavelet.hpp, LaplaceUpGradientPrewavelet.hpp and
 * LaplaceDownPrewavelet.hpp.
 */
class OperationLaplacePrewavelet: public UpDownOneOpDimWithShadow
{
public:
	/**
	 * Construtor of OperationLaplacePrewavelet
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationLaplacePrewavelet(GridStorage* storage, GridStorage* shadowstorage) :
			UpDownOneOpDimWithShadow(storage,shadowstorage)
	{

	}

	/**
	 * Destructor
	 */
	virtual ~OperationLaplacePrewavelet()
	{
	}


protected:

	virtual void up(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::LaplaceUpPrewavelet func(this->storage);
		sweep<detail::LaplaceUpPrewavelet> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

	virtual void down(DataVector& alpha, DataVector& result, size_t dim)
	{
	}

	virtual void downOpDim(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::LaplaceDownGradientPrewavelet func(this->storage);
		sweep<detail::LaplaceDownGradientPrewavelet> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

	virtual void upOpDim(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::LaplaceUpGradientPrewavelet func(this->storage);
		sweep<detail::LaplaceUpGradientPrewavelet> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

};

}

#endif /* OPERATIONLAPLACEPREWAVELET_HPP */
