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
