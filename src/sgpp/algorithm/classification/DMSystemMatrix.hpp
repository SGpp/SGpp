/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#ifndef DMSYSTEMMATRIX_HPP
#define DMSYSTEMMATRIX_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/OperationB.hpp"
#include "operation/OperationMatrix.hpp"

namespace sg
{

/**
 * Class that implements the virtual class OperationMatrix for the
 * application of classification for the Systemmatrix
 */
class DMSystemMatrix : public OperationMatrix
{
private:
	/// the lambda, the regularisation parameter
	double lamb;
	/// OperationMatrix, the regularisation mehtod
	OperationMatrix* C;
	/// OperationB for calculating the data matrix
	OperationB* B;
	/// Pointer to the data vector
	DataVector* data;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to DataVector that contains the training data
	 * @param C the regression functional
	 * @param lambda the lambda, the regression parameter
	 */
	DMSystemMatrix(Grid& SparseGrid, DataVector& trainData, OperationMatrix& C, double lambda);

	/**
	 * Std-Destructor
	 */
	virtual ~DMSystemMatrix();

	virtual void mult(DataVector& alpha, DataVector& result);

	/**
	 * Generates the right hand side of the classification equation
	 *
	 * @param classes the class information of the training data
	 * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
	 */
	void generateb(DataVector& classes, DataVector& b);
};

}

#endif /* DMSYSTEMMATRIX_HPP */
