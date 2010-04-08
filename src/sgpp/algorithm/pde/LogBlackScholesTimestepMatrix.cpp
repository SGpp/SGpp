/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*               2010      Stefanie Schraufstetter (schraufs@in.tum.de)      */
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

#include "algorithm/pde/LogBlackScholesTimestepMatrix.hpp"
#include "exception/algorithm_exception.hpp"
#include <cmath>

namespace sg
{

LogBlackScholesTimestepMatrix::LogBlackScholesTimestepMatrix(Grid& SparseGrid, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode)
{
	this->myGrid = &SparseGrid;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->r = r;
	this->mus = &mu;
	this->sigmas = &sigma;
	this->rhos = &rho;

	// build the coefficient vectors for the operations
	this->gammaCoef = new DataVector(SparseGrid.getStorage()->dim(), SparseGrid.getStorage()->dim());
	this->deltaCoef = new DataVector(SparseGrid.getStorage()->dim());
	buildDeltaCoefficients();
	buildGammaCoefficients();

	// Create needed operations
	this->OpDelta = SparseGrid.createOperationDeltaLog(*this->deltaCoef);
	this->OpGamma = SparseGrid.createOperationGammaLog(*this->gammaCoef);
	this->OpLTwo = SparseGrid.createOperationLTwoDotProduct();
}

LogBlackScholesTimestepMatrix::~LogBlackScholesTimestepMatrix()
{
	delete this->OpDeltaLog;
	delete this->OpGammaLog;
	delete this->OpLTwo;
	delete this->gammaCoef;
	delete this->deltaCoef;
	delete this->BoundaryUpdate;
}

void LogBlackScholesTimestepMatrix::mult(DataVector& alpha, DataVector& result)
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

		applyMassMatrix(alpha, temp);
		result.add_parallel(temp);

		applyLOperator(alpha, temp);
		result.axpy_parallel((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrix(alpha, temp);
		result.add_parallel(temp);

		applyLOperator(alpha, temp);
		result.axpy_parallel((-0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("LogBlackScholesTimestepMatrix::mult : An unknown operation mode was specified!");
	}
}

void LogBlackScholesTimestepMatrix::generateRHS(DataVector& data, DataVector& rhs)
{
	if (this->tOperationMode == "ExEul")
	{
		DataVector temp(data.getSize());

		applyMassMatrix(data, temp);
		rhs.add_parallel(temp);

		applyLOperator(data, temp);
		rhs.axpy_parallel(this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs.setAll(0.0);

		applyMassMatrix(data, rhs);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs.setAll(0.0);

		DataVector temp(data.getSize());

		applyMassMatrix(data, temp);
		rhs.add_parallel(temp);

		applyLOperator(data, temp);
		rhs.axpy_parallel((0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("LogBlackScholesTimestepMatrix::generateRHS : An unknown operation mode was specified!");
	}
}

void LogBlackScholesTimestepMatrix::applyLOperator(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpLTwo->mult(alpha, temp);
		result.axpy_parallel((-1.0)*this->r, temp);
	}

	// Apply the delta method
	this->OpDeltaLog->mult(alpha, temp);
	result.add_parallel(temp);

	// Apply the gamma method
	this->OpGammaLog->mult(alpha, temp);
	result.sub_parallel(temp);
}

void LogBlackScholesTimestepMatrix::applyMassMatrix(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwo->mult(alpha, temp);

	result.add_parallel(temp);
}

void LogBlackScholesTimestepMatrix::finishTimestep(DataVector& alpha)
{
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "ExEul")
		{
			this->BoundaryUpdate->multiplyBoundary(alpha, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
}

void LogBlackScholesTimestepMatrix::startTimestep(DataVector& alpha)
{
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{
			this->BoundaryUpdate->multiplyBoundary(alpha, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
}

Grid* LogBlackScholesTimestepMatrix::getGrid()
{
	return myGrid;
}

void LogBlackScholesTimestepMatrix::buildGammaCoefficients()
{
	size_t dim = this->myGrid->getStorage()->dim();

	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
				this->gammaCoef->set((dim*i)+j, 0.5*((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get((i*dim)+j)));
			}
			else
			{
				this->gammaCoef->set((dim*i)+j, ((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get((i*dim)+j)));
			}
		}
	}
}

void LogBlackScholesTimestepMatrix::buildDeltaCoefficients()
{
	size_t dim = this->myGrid->getStorage()->dim();
	double covar_sum = 0.0;

	for (size_t i = 0; i < dim; i++)
	{
		covar_sum = 0.0;
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
				covar_sum += ((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get((i*dim)+j));
			}
			else
			{
				covar_sum += (0.5*((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get((i*dim)+j)));
			}
		}
		this->deltaCoef->set(i, this->mus->get(i)-covar_sum);
	}
}

}
