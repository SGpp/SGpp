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
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;
	this->InnerGrid = NULL;
	this->alpha_inner = NULL;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->GridConverter = new DirichletGridConverter();
	this->r = r;
	this->mus = &mu;
	this->sigmas = &sigma;
	this->rhos = &rho;

	// build the coefficient vectors for the operations
	this->gammaCoef = new DataVector(SparseGrid.getStorage()->dim(), SparseGrid.getStorage()->dim());
	this->deltaCoef = new DataVector(SparseGrid.getStorage()->dim());
	buildDeltaCoefficients();
	buildGammaCoefficients();

	// Create needed operations, on boundary grid
	this->OpDeltaBound = this->BoundGrid->createOperationDelta(*this->deltaCoef);
	this->OpGammaBound = this->BoundGrid->createOperationGamma(*this->gammaCoef);
	this->OpLTwoBound = this->BoundGrid->createOperationLTwoDotProduct();

	// create the inner grid
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

	//Create needed operations, on inner grid
	this->OpDeltaInner = this->InnerGrid->createOperationDelta(*this->deltaCoef);
	this->OpGammaInner = this->InnerGrid->createOperationGamma(*this->gammaCoef);
	this->OpLTwoInner = this->InnerGrid->createOperationLTwoDotProduct();

	// right hand side if System
	this->rhs = new DataVector(1);
}

BlackScholesODESolverSystem::~BlackScholesODESolverSystem()
{
	delete this->OpDeltaBound;
	delete this->OpGammaBound;
	delete this->OpLTwoBound;
	delete this->OpDeltaInner;
	delete this->OpGammaInner;
	delete this->OpLTwoInner;
	delete this->gammaCoef;
	delete this->deltaCoef;
	delete this->BoundaryUpdate;
	delete this->GridConverter;
	if (this->InnerGrid != NULL)
	{
		delete this->InnerGrid;
	}
	if (this->alpha_inner != NULL)
	{
		delete this->alpha_inner;
	}
	delete this->rhs;
}

void BlackScholesODESolverSystem::mult(DataVector& alpha, DataVector& result)
{
//	std::cout << "entered mult" << std::endl;
//	std::cout << result.getSize() << std::endl;
//	std::cout << alpha.getSize() << std::endl;

	if (this->tOperationMode == "ExEul")
	{
		result.setAll(0.0);

		applyMassMatrixInner(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrixInner(alpha, temp);
		result.add_parallel(temp);

		applyLOperatorInner(alpha, temp);
		result.axpy_parallel((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrixInner(alpha, temp);
		result.add_parallel(temp);

		applyLOperatorInner(alpha, temp);
		result.axpy_parallel((-0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("BlackScholesTimestepMatrix::mult : An unknown operation mode was specified!");
	}
}

DataVector* BlackScholesODESolverSystem::generateRHS()
{
	//std::cout << "Entered generateRHS" << std::endl;

	DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add_parallel(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy_parallel(this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs_complete.setAll(0.0);

		applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add_parallel(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy_parallel((0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("BlackScholesTimestepMatrix::generateRHS : An unknown operation mode was specified!");
	}

	// Now we have the right hand side, lets apply the riskfree rate for the next timestep
	this->startTimestep();

	// Now apply the boundary ansatzfunctions to the inner ansatzfunctions
	DataVector result_complete(this->alpha_complete->getSize());
	DataVector alpha_bound(*this->alpha_complete);

	result_complete.setAll(0.0);

	this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

	// apply CG Matrix
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else if (this->tOperationMode == "ImEul")
	{
		DataVector temp(alpha_bound.getSize());
		temp.setAll(0.0);

		applyMassMatrixComplete(alpha_bound, temp);
		result_complete.add_parallel(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy_parallel((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		DataVector temp(alpha_bound.getSize());
		temp.setAll(0.0);

		applyMassMatrixComplete(alpha_bound, temp);
		result_complete.add_parallel(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy_parallel((-0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("BlackScholesTimestepMatrix::generateRHS : An unknown operation mode was specified (mult)!");
	}
	rhs_complete.sub(result_complete);

	this->rhs->resize(this->alpha_inner->getSize());

	//std::cout << this->rhs->getSize() << std::endl;
	//std::cout << rhs_complete.getSize() << std::endl;

	this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

	return this->rhs;
}

void BlackScholesODESolverSystem::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpLTwoBound->mult(alpha, temp);
		result.axpy_parallel((-1.0)*this->r, temp);
	}

	// Apply the delta method
	this->OpDeltaBound->mult(alpha, temp);
	result.add_parallel(temp);

	// Apply the gamma method
	this->OpGammaBound->mult(alpha, temp);
	result.sub_parallel(temp);
}

void BlackScholesODESolverSystem::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpLTwoInner->mult(alpha, temp);
		result.axpy_parallel((-1.0)*this->r, temp);
	}

	// Apply the delta method
	this->OpDeltaInner->mult(alpha, temp);
	result.add_parallel(temp);

	// Apply the gamma method
	this->OpGammaInner->mult(alpha, temp);
	result.sub_parallel(temp);
}

void BlackScholesODESolverSystem::applyMassMatrixComplete(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add_parallel(temp);
}

void BlackScholesODESolverSystem::applyMassMatrixInner(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoInner->mult(alpha, temp);

	result.add_parallel(temp);
}

void BlackScholesODESolverSystem::finishTimestep()
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "ExEul")
		{
			this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
}

void BlackScholesODESolverSystem::startTimestep()
{
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{
			this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
}

Grid* BlackScholesODESolverSystem::getGrid()
{
	return BoundGrid;
}

void BlackScholesODESolverSystem::buildGammaCoefficients()
{
	size_t dim = this->BoundGrid->getStorage()->dim();

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
	size_t dim = this->BoundGrid->getStorage()->dim();
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

DataVector* BlackScholesODESolverSystem::getGridCoefficientsForCG()
{
	this->GridConverter->calcInnerCoefs(*this->alpha_complete, *this->alpha_inner);

	return this->alpha_inner;
}

DataVector* BlackScholesODESolverSystem::getGridCoefficients()
{
	return this->alpha_complete;
}

}
