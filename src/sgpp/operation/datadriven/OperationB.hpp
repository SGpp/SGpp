/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008-2009 Dirk Pflueger (pflueged@in.tum.de)                */
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

#ifndef OPERATIONB_HPP
#define OPERATIONB_HPP

#include "data/DataVector.hpp"

namespace sg
{

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 * 
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points, then B is a @f$N\times m@f$ matrix, with
 * @f[ (B)_{i,j} = \varphi_i(x_j). @f]
 */
class OperationB
{
public:
	/**
	 * Constructor
	 */
	OperationB() {}

	/**
	 * Destructor
	 */
	virtual ~OperationB() {}

	/**
	 * Multiplication of @f$B@f$ with vector @f$\alpha@f$
	 *
	 * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
	 * @param data vector, providing the data points x row-wise
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void mult(DataVector& alpha, DataVector& data, DataVector& result) = 0;

	/**
	 * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
	 *
	 * @param alpha vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
	 * @param data vector, providing the data points x row-wise
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void multTranspose(DataVector& alpha, DataVector& data, DataVector& result) = 0;
};

}

#endif /* OPERATIONB_HPP */
