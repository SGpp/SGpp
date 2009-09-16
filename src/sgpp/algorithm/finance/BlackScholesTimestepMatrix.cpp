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

#include "algorithm/finance/BlackScholesTimestepMatrix.hpp"

namespace sg
{

BlackScholesTimestepMatrix::BlackScholesTimestepMatrix(Grid& SparseGrid, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, bool bCrankNicolsonMatrix)
{
	this->OpDelta = SparseGrid.createOperationDelta(mu);
	this->OpGammaOne = SparseGrid.createOperationGammaPartOne(sigma, rho);
	this->OpGammaTwo = SparseGrid.createOperationGammaPartTwo(sigma, rho);
	this->OpGammaThree = SparseGrid.createOperationGammaPartThree(sigma, rho);
	this->OpRiskfree = SparseGrid.createOperationRiskfreeRate();
	this->r = r;
	this->mus = &mu;
	this->sigmas = &sigma;
	this->rhos = &rho;
	this->bIsCrankNicolsonMatrix = bCrankNicolsonMatrix;
	this->TimestepSize = TimestepSize;
	this->myGrid = &SparseGrid;
}

BlackScholesTimestepMatrix::~BlackScholesTimestepMatrix()
{
	delete this->OpDelta;
	delete this->OpGammaOne;
	delete this->OpGammaTwo;
	delete this->OpGammaThree;
	delete this->OpRiskfree;
}

void BlackScholesTimestepMatrix::mult(DataVector& alpha, DataVector& result)
{
	if (this->bIsCrankNicolsonMatrix == false)
	{
		applyL(alpha, result);
	}
	else
	{
		result.setAll(0.0);

		applyL(alpha, result);
		result.mult((-0.5)*this->TimestepSize);
		result.add(alpha);
	}
}

void BlackScholesTimestepMatrix::generateRHS(DataVector& data, DataVector& rhs)
{
	if (this->bIsCrankNicolsonMatrix == false)
	{
		// @todo (heinecke) throw some exception here
	}
	else
	{
		rhs.setAll(0.0);

		applyL(data, rhs);

		rhs.mult(0.5*this->TimestepSize);
		rhs.add(data);
	}
}

void BlackScholesTimestepMatrix::applyL(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	// Apply the riskfree rate
	this->OpRiskfree->mult(alpha, temp);
	result.axpy((-1.0)*this->r, temp);

	// Apply the delta method
	this->OpDelta->mult(alpha, temp);
	result.add(temp);

	// Apply the gamma method, part 3
	this->OpGammaThree->mult(alpha, temp);
	result.sub(temp);

	// Apply the gamma method, part 2
	this->OpGammaTwo->mult(alpha, temp);
	result.sub(temp);

	// Apply the gamma method, part 1
	this->OpGammaOne->mult(alpha, temp);
	result.add(temp);
}

Grid* BlackScholesTimestepMatrix::getGrid()
{
	return myGrid;
}

}
