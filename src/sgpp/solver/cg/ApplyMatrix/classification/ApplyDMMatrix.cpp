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

#include "solver/cg/ApplyMatrix/classification/ApplyDMMatrix.hpp"
#include "exception/operation_exception.hpp"


namespace sg
{

ApplyDMMatrix::ApplyDMMatrix(Grid* SparseGrid, std::string StiffnessMode, double lambda)
{
	this->StiffMode = StiffnessMode;

	if (this->StiffMode != "L" && this->StiffMode != "I")
	{
		throw new operation_exception("You have chosen an invalid stiffness mode in ApplyDMMatrix! L or I are valid");
	}
	// create the operations needed in ApplyMatrix
	this->C = SparseGrid->createOperationLaplace();
	this->B = SparseGrid->createOperationB();
	this->lamb = lambda;
}

ApplyDMMatrix::~ApplyDMMatrix()
{
	delete this->C;
	delete this->B;
}

void ApplyDMMatrix::operator()(DataVector& data, DataVector& alpha, DataVector& result)
{
	DataVector temp(data.getSize());
    size_t M = data.getSize();

    // Operation B
    this->B->multTranspose(alpha, data, temp);
    this->B->mult(temp, data, result);

    // Laplace Matrix
    if (this->StiffMode == "L")
    {
		DataVector temptwo(alpha.getSize());
		this->C->mult(alpha, temptwo);
		result.axpy(M*this->lamb, temptwo);
    }

    // Identity Matrix
    if (this->StiffMode == "I")
    {
		result.axpy(M*this->lamb, alpha);
    }
}

void ApplyDMMatrix::generateb(DataVector& data, DataVector& classes, DataVector& b)
{
	this->B->mult(classes, data, b);
}

}
