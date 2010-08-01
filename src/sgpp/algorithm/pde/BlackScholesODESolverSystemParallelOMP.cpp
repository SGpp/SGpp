/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/BlackScholesODESolverSystemParallelOMP.hpp"
#include "exception/algorithm_exception.hpp"

#include "algorithm/pde/StdUpDown.hpp"
#include "algorithm/pde/UpDownOneOpDim.hpp"
#include "algorithm/pde/UpDownTwoOpDims.hpp"

#ifdef USEOMPTHREE
#include "omp.h"
#endif

namespace sg
{

BlackScholesODESolverSystemParallelOMP::BlackScholesODESolverSystemParallelOMP(Grid& SparseGrid, DataVector& alpha, DataVector& mu,
			DataVector& sigma, DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
			bool bLogTransform, bool useCoarsen, double coarsenThreshold, double coarsenPercent,
			size_t numExecCoarsen, size_t MPIRank) : BlackScholesODESolverSystem(SparseGrid, alpha, mu, sigma, rho,
			r, TimestepSize, OperationMode, bLogTransform, useCoarsen, coarsenThreshold, coarsenPercent, numExecCoarsen, MPIRank)
{}

BlackScholesODESolverSystemParallelOMP::~BlackScholesODESolverSystemParallelOMP()
{
}

void BlackScholesODESolverSystemParallelOMP::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

#ifdef USEOMPTHREE
	size_t nDims = this->InnerGrid->getStorage()->dim();
	omp_lock_t resultMutex;
	omp_init_lock(&resultMutex);

	#pragma omp parallel shared(resultMutex, alpha, result)
	{
#ifndef AIX_XLC
		#pragma omp single nowait
#endif
#ifdef AIX_XLC
		#pragma omp single
#endif
		{
			// Apply the riskfree rate
			#pragma omp task shared(alpha, result, resultMutex)
			{
				if (this->r != 0.0)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((StdUpDown*)(this->OpLTwoInner))->multParallelBuildingBlock(alpha, myResult);

					// semaphore
					omp_set_lock(&resultMutex);
					result.axpy((-1.0)*this->r, myResult);
					omp_unset_lock(&resultMutex);
				}
			}

			// Apply the delta method
			for (size_t i = 0; i < nDims; i++)
			{
				#pragma omp task firstprivate(i) shared(alpha, result, resultMutex)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownOneOpDim*)(this->OpDeltaInner))->multParallelBuildingBlock(alpha, myResult, i);

					// semaphore
					omp_set_lock(&resultMutex);
					result.add(myResult);
					omp_unset_lock(&resultMutex);
				}
			}

			// Apply the gamma method
			for (size_t i = 0; i < nDims; i++)
			{
				for (size_t j = 0; j < nDims; j++)
				{
					// symmetric
					if (j <= i)
					{
						#pragma omp task firstprivate(i, j) shared(alpha, result, resultMutex)
						{
							DataVector myResult(result.getSize());

							/// @todo (heinecke) discuss methods in order to avoid this cast
							((UpDownTwoOpDims*)(this->OpGammaInner))->multParallelBuildingBlock(alpha, myResult, i, j);

							// semaphore
							omp_set_lock(&resultMutex);
							result.sub(myResult);
							omp_unset_lock(&resultMutex);
						}
					}
				}
			}

			#pragma omp taskwait
		}
	}

	omp_destroy_lock(&resultMutex);
#else
	DataVector temp(alpha.getSize());

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpLTwoInner->mult(alpha, temp);
		result.axpy((-1.0)*this->r, temp);
	}

	// Apply the delta method
	this->OpDeltaInner->mult(alpha, temp);
	result.add(temp);

	// Apply the gamma method
	this->OpGammaInner->mult(alpha, temp);
	result.sub(temp);
#endif
}

void BlackScholesODESolverSystemParallelOMP::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

#ifdef USEOMPTHREE
	size_t nDims = this->InnerGrid->getStorage()->dim();
	omp_lock_t resultMutex;
	omp_init_lock(&resultMutex);

	#pragma omp parallel shared(resultMutex, alpha, result)
	{
#ifndef AIX_XLC
		#pragma omp single nowait
#endif
#ifdef AIX_XLC
		#pragma omp single
#endif
		{
			// Apply the riskfree rate
			#pragma omp task shared(alpha, result, resultMutex)
			{
				if (this->r != 0.0)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((StdUpDown*)(this->OpLTwoBound))->multParallelBuildingBlock(alpha, myResult);

					// semaphore
					omp_set_lock(&resultMutex);
					result.axpy((-1.0)*this->r, myResult);
					omp_unset_lock(&resultMutex);
				}
			}

			// Apply the delta method
			for (size_t i = 0; i < nDims; i++)
			{
				#pragma omp task firstprivate(i) shared(alpha, result, resultMutex)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownOneOpDim*)(this->OpDeltaBound))->multParallelBuildingBlock(alpha, myResult, i);

					// semaphore
					omp_set_lock(&resultMutex);
					result.add(myResult);
					omp_unset_lock(&resultMutex);
				}
			}

			// Apply the gamma method
			for (size_t i = 0; i < nDims; i++)
			{
				for (size_t j = 0; j < nDims; j++)
				{
					// symmetric
					if (j <= i)
					{
						#pragma omp task firstprivate(i, j) shared(alpha, result, resultMutex)
						{
							DataVector myResult(result.getSize());

							/// @todo (heinecke) discuss methods in order to avoid this cast
							((UpDownTwoOpDims*)(this->OpGammaBound))->multParallelBuildingBlock(alpha, myResult, i, j);

							// semaphore
							omp_set_lock(&resultMutex);
							result.sub(myResult);
							omp_unset_lock(&resultMutex);
						}
					}
				}
			}

			#pragma omp taskwait
		}
	}

	omp_destroy_lock(&resultMutex);
#else
	DataVector temp(alpha.getSize());

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpLTwoBound->mult(alpha, temp);
		result.axpy((-1.0)*this->r, temp);
	}

	// Apply the delta method
	this->OpDeltaBound->mult(alpha, temp);
	result.add(temp);

	// Apply the gamma method
	this->OpGammaBound->mult(alpha, temp);
	result.sub(temp);
#endif
}

void BlackScholesODESolverSystemParallelOMP::applyMassMatrixInner(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoInner->mult(alpha, temp);

	result.add(temp);
}

void BlackScholesODESolverSystemParallelOMP::applyMassMatrixComplete(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

}
