/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
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

#ifndef OPERATIONLAPLACELINEAR_HPP
#define OPERATIONLAPLACELINEAR_HPP

#include "operation/common/OperationMatrix.hpp"

#include "algorithm/datadriven/UnidirGradient.hpp"

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg
{

/**
 * Implementation for linear functions of Laplace Operation, linear grids without boundaries
 *
 * @version $HEAD$
 */
class OperationLaplaceLinear: public OperationMatrix, public UnidirGradient
{
public:
	/**
	 * Construtor of OperationLaplaceLinear
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationLaplaceLinear(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~OperationLaplaceLinear();

	virtual void mult(DataVector& alpha, DataVector& result);

protected:
#ifndef USEOMPTHREE
	virtual void gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);
#endif

#ifdef USEOMPTHREE
	virtual void gradient_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);
#endif

	virtual void up(DataVector& alpha, DataVector& result, size_t dim);

	virtual void down(DataVector& alpha, DataVector& result, size_t dim);

	virtual void downGradient(DataVector& alpha, DataVector& result, size_t dim);

	virtual void upGradient(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* OPERATIONLAPLACELINEAR_HPP */
