/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp"
#include "exception/algorithm_exception.hpp"

#include "algorithm/pde/StdUpDown.hpp"
#include "algorithm/pde/UpDownOneOpDim.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace sg
{

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs)
{
	this->Laplace_Complete = this->BoundGrid->createOperationLaplace();
	this->Laplace_Inner = this->InnerGrid->createOperationLaplace();
}

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI()
{
	delete this->Laplace_Complete;
	delete this->Laplace_Inner;
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	#pragma omp parallel shared(alpha, result)
	{
		#pragma omp single nowait
		{
			std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
			size_t nDims = algoDims.size();
#ifdef _OPENMP
			omp_lock_t Mutex;
			omp_init_lock(&Mutex);
#endif
			// Apply Laplace, parallel in Dimensions
			for (size_t i = 0; i < nDims; i++)
			{
				#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownOneOpDim*)(this->Laplace_Inner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

					// semaphore
#ifdef _OPENMP
					omp_set_lock(&Mutex);
#endif
					result.add(myResult);
#ifdef _OPENMP
					omp_unset_lock(&Mutex);
#endif
				}
			}

			#pragma omp taskwait

#ifdef _OPENMP
			omp_destroy_lock(&Mutex);
#endif
		}
	}
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	#pragma omp parallel shared(alpha, result)
	{
		#pragma omp single nowait
		{
			std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
			size_t nDims = algoDims.size();
#ifdef _OPENMP
			omp_lock_t Mutex;
			omp_init_lock(&Mutex);
#endif
			// Apply Laplace, parallel in Dimensions
			for (size_t i = 0; i < nDims; i++)
			{
				#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownOneOpDim*)(this->Laplace_Complete))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

					// semaphore
#ifdef _OPENMP
					omp_set_lock(&Mutex);
#endif
					result.add(myResult);
#ifdef _OPENMP
					omp_unset_lock(&Mutex);
#endif
				}
			}

			#pragma omp taskwait

#ifdef _OPENMP
			omp_destroy_lock(&Mutex);
#endif
		}
	}
}

}
