/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "algorithm/pde/LogBlackScholesODESolverSystem.hpp"
#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include <cmath>

namespace sg
{

LogBlackScholesODESolverSystem::LogBlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode)
	: BlackScholesODESolverSystem(SparseGrid, alpha, mu, sigma, rho, r, TimestepSize, OperationMode)
{
	// create needed operations that are different in case of a log-transformed Black-Scholoes equation

	// operations on boundary grid
	this->OpDeltaBound = this->BoundGrid->createOperationDeltaLog(*this->deltaCoef);
	this->OpGammaBound = this->BoundGrid->createOperationGammaLog(*this->gammaCoef);

	//operations on inner grid
	this->OpDeltaInner = this->InnerGrid->createOperationDeltaLog(*this->deltaCoef);
	this->OpGammaInner = this->InnerGrid->createOperationGammaLog(*this->gammaCoef);

}




void LogBlackScholesODESolverSystem::buildGammaCoefficients()
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

void LogBlackScholesODESolverSystem::buildDeltaCoefficients()
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

}
