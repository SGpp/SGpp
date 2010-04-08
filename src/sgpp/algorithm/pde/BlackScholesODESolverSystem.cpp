/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include <cmath>

namespace sg
{

BlackScholesODESolverSystem::BlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode)
{
	this->myGrid = &SparseGrid;
	this->alpha_complete = &alpha;
	this->alpha_inner = &alpha;
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
	this->OpDelta = SparseGrid.createOperationDelta(*this->deltaCoef);
	this->OpGamma = SparseGrid.createOperationGamma(*this->gammaCoef);
	this->OpLTwo = SparseGrid.createOperationLTwoDotProduct();
}

BlackScholesODESolverSystem::~BlackScholesODESolverSystem()
{
	delete this->OpDelta;
	delete this->OpGamma;
	delete this->OpLTwo;
	delete this->gammaCoef;
	delete this->deltaCoef;
	delete this->BoundaryUpdate;
}

void BlackScholesODESolverSystem::mult(DataVector& alpha, DataVector& result)
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
		throw new algorithm_exception("BlackScholesTimestepMatrix::mult : An unknown operation mode was specified!");
	}
}

void BlackScholesODESolverSystem::generateRHS(DataVector& data, DataVector& rhs)
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
		throw new algorithm_exception("BlackScholesTimestepMatrix::generateRHS : An unknown operation mode was specified!");
	}
}

void BlackScholesODESolverSystem::applyLOperator(DataVector& alpha, DataVector& result)
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
	this->OpDelta->mult(alpha, temp);
	result.add_parallel(temp);

	// Apply the gamma method
	this->OpGamma->mult(alpha, temp);
	result.sub_parallel(temp);
}

void BlackScholesODESolverSystem::applyMassMatrix(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwo->mult(alpha, temp);

	result.add_parallel(temp);
}

void BlackScholesODESolverSystem::finishTimestep(DataVector& alpha)
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

void BlackScholesODESolverSystem::startTimestep(DataVector& alpha)
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

Grid* BlackScholesODESolverSystem::getGrid()
{
	return myGrid;
}

void BlackScholesODESolverSystem::buildGammaCoefficients()
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

void BlackScholesODESolverSystem::buildDeltaCoefficients()
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

DataVector* BlackScholesODESolverSystem::getGridCoefficients()
{
	return this->alpha_complete;
}

DataVector* BlackScholesODESolverSystem::getGridCoefficientsForCG()
{
	return this->alpha_inner;
}

}
