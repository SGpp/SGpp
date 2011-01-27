/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichletParallelOMP.hpp"
#include "exception/algorithm_exception.hpp"

#include "algorithm/pde/StdUpDown.hpp"
#include "algorithm/pde/UpDownOneOpDim.hpp"

namespace sg
{

PoissonEquationEllipticPDESolverSystemDirichletParallelOMP::PoissonEquationEllipticPDESolverSystemDirichletParallelOMP(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs)
{
	this->Laplace_Complete = this->BoundGrid->createOperationLaplace();
	this->Laplace_Inner = this->InnerGrid->createOperationLaplace();
}

PoissonEquationEllipticPDESolverSystemDirichletParallelOMP::~PoissonEquationEllipticPDESolverSystemDirichletParallelOMP()
{
	delete this->Laplace_Complete;
	delete this->Laplace_Inner;
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelOMP::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

#ifdef USEOMPTHREE
	#pragma omp parallel shared(alpha, result)
	{
	#ifndef AIX_XLC
		#pragma omp single nowait
	#endif
	#ifdef AIX_XLC
		#pragma omp single
	#endif
		{
			std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
			size_t nDims = algoDims.size();
			omp_lock_t Mutex;
			omp_init_lock(&Mutex);

			// Apply Laplace, parallel in Dimensions
			for (size_t i = 0; i < nDims; i++)
			{
				#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownOneOpDim*)(this->Laplace_Inner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

					// semaphore
					omp_set_lock(&Mutex);
					result.add(myResult);
					omp_unset_lock(&Mutex);
				}
			}

			#pragma omp taskwait

			omp_destroy_lock(&Mutex);
		}
	}
#else
	Laplace_Inner->mult(alpha, result);
#endif
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelOMP::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

#ifdef USEOMPTHREE
	#pragma omp parallel shared(alpha, result)
	{
	#ifndef AIX_XLC
		#pragma omp single nowait
	#endif
	#ifdef AIX_XLC
		#pragma omp single
	#endif
		{
			std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
			size_t nDims = algoDims.size();
			omp_lock_t Mutex;
			omp_init_lock(&Mutex);

			// Apply Laplace, parallel in Dimensions
			for (size_t i = 0; i < nDims; i++)
			{
				#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownOneOpDim*)(this->Laplace_Complete))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

					// semaphore
					omp_set_lock(&Mutex);
					result.add(myResult);
					omp_unset_lock(&Mutex);
				}
			}

			#pragma omp taskwait

			omp_destroy_lock(&Mutex);
		}
	}
#else
	Laplace_Complete->mult(alpha, result);
#endif
}

}
