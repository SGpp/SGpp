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

#include "algorithm/pde/BlackScholesTimestepMatrix.hpp"

namespace sg
{

BlackScholesTimestepMatrix::BlackScholesTimestepMatrix(Grid& SparseGrid, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode)
{
	this->OpDelta = SparseGrid.createOperationDelta(mu);
	this->OpGammaOne = SparseGrid.createOperationGammaPartOne(sigma, rho);
	this->OpGammaTwo = SparseGrid.createOperationGammaPartTwo(sigma, rho);
	this->OpGammaThree = SparseGrid.createOperationGammaPartThree(sigma, rho);
	this->OpRiskfree = SparseGrid.createOperationRiskfreeRate();
	this->OpMass = SparseGrid.createOperationLTwoDotProduct();
	this->r = r;
	this->mus = &mu;
	this->sigmas = &sigma;
	this->rhos = &rho;
	this->tOperationMode = OperationMode;
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
	delete this->OpMass;
}

void BlackScholesTimestepMatrix::mult(DataVector& alpha, DataVector& result)
{
	if (this->tOperationMode == "ExEul")
	{
		result.setAll(0.0);

		applyMassMatrix(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());
		temp.setAll(0.0);

		applyMassMatrix(alpha, temp);
		result.add(temp);

		temp.setAll(0.0);
		applyLOperator(alpha, temp);
		result.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		/*result.setAll(0.0);

		applyL(alpha, result);
		result.mult((-0.5)*this->TimestepSize);
		result.add(alpha);*/
	}
	else
	{
		// @todo (heinecke) throw some exception here
	}
}

void BlackScholesTimestepMatrix::generateRHS(DataVector& data, DataVector& rhs)
{
	if (this->tOperationMode == "ExEul")
	{
		DataVector temp(data.getSize());
		temp.setAll(0.0);

		applyMassMatrix(data, temp);
		rhs.add(temp);

		temp.setAll(0.0);
		applyLOperator(data, temp);
		rhs.axpy(this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "ImEul")
	{
		applyMassMatrix(data, rhs);
	}
	else if (this->tOperationMode == "CrNic")
	{
		/*rhs.setAll(0.0);

		applyL(data, rhs);

		rhs.mult(0.5*this->TimestepSize);
		rhs.add(data);*/
	}
	else
	{
		// @todo (heinecke) throw some exception here
	}
}

void BlackScholesTimestepMatrix::applyLOperator(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpRiskfree->mult(alpha, temp);
		result.axpy((-1.0)*this->r, temp);
	}

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
	//this->OpGammaOne->mult(alpha, temp);
	//result.add(temp);
}

void BlackScholesTimestepMatrix::applyMassMatrix(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	// Apply the mass matrix
	this->OpMass->mult(alpha, temp);

	result.add(temp);
}

Grid* BlackScholesTimestepMatrix::getGrid()
{
	return myGrid;
}

}
