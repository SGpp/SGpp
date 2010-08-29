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
			bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel) : BlackScholesODESolverSystem(SparseGrid, alpha, mu, sigma, rho,
			r, TimestepSize, OperationMode, bLogTransform, useCoarsen, coarsenThreshold, adaptSolveMode, numCoarsenPoints, refineThreshold, refineMode, refineMaxLevel)
{}

BlackScholesODESolverSystemParallelOMP::~BlackScholesODESolverSystemParallelOMP()
{
}

void BlackScholesODESolverSystemParallelOMP::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

#ifdef USEOMPTHREE
	std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();
	omp_lock_t DeltaMutex;
	omp_lock_t GammaMutex;
	omp_init_lock(&DeltaMutex);
	omp_init_lock(&GammaMutex);
	DataVector DeltaResult(result);
	DataVector GammaResult(result);

	// Apply the riskfree rate
	#pragma omp task shared(alpha, result)
	{
		if (this->r != 0.0)
		{
			DataVector myResult(result.getSize());

			/// @todo (heinecke) discuss methods in order to avoid this cast
			((StdUpDown*)(this->OpLTwoInner))->multParallelBuildingBlock(alpha, myResult);

			// no semaphore needed
			result.axpy((-1.0)*this->r, myResult);
		}
	}

	// Apply the delta method
	for (size_t i = 0; i < nDims; i++)
	{
		#pragma omp task firstprivate(i) shared(alpha, DeltaMutex, DeltaResult, result)
		{
			DataVector myResult(result.getSize());

			/// @todo (heinecke) discuss methods in order to avoid this cast
			((UpDownOneOpDim*)(this->OpDeltaInner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

			// semaphore
			omp_set_lock(&DeltaMutex);
			DeltaResult.add(myResult);
			omp_unset_lock(&DeltaMutex);
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
				#pragma omp task firstprivate(i, j) shared(alpha, GammaMutex, GammaResult, result)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownTwoOpDims*)(this->OpGammaInner))->multParallelBuildingBlock(alpha, myResult, algoDims[i], algoDims[j]);

					// semaphore
					omp_set_lock(&GammaMutex);
					GammaResult.add(myResult);
					omp_unset_lock(&GammaMutex);
				}
			}
		}
	}

	#pragma omp taskwait

	omp_destroy_lock(&GammaMutex);
	omp_destroy_lock(&DeltaMutex);

	// sum up
	result.add(DeltaResult);
	result.sub(GammaResult);
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
	std::vector<size_t> algoDims = this->BoundGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();
	omp_lock_t DeltaMutex;
	omp_lock_t GammaMutex;
	omp_init_lock(&DeltaMutex);
	omp_init_lock(&GammaMutex);
	DataVector DeltaResult(result);
	DataVector GammaResult(result);

	// Apply the riskfree rate
	#pragma omp task shared(alpha, result)
	{
		if (this->r != 0.0)
		{
			DataVector myResult(result.getSize());

			/// @todo (heinecke) discuss methods in order to avoid this cast
			((StdUpDown*)(this->OpLTwoBound))->multParallelBuildingBlock(alpha, myResult);

			// no semaphore needed
			result.axpy((-1.0)*this->r, myResult);
		}
	}

	// Apply the delta method
	for (size_t i = 0; i < nDims; i++)
	{
		#pragma omp task firstprivate(i) shared(alpha, DeltaMutex, DeltaResult, result)
		{
			DataVector myResult(result.getSize());

			/// @todo (heinecke) discuss methods in order to avoid this cast
			((UpDownOneOpDim*)(this->OpDeltaBound))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

			// semaphore
			omp_set_lock(&DeltaMutex);
			DeltaResult.add(myResult);
			omp_unset_lock(&DeltaMutex);
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
				#pragma omp task firstprivate(i, j) shared(alpha, GammaMutex, GammaResult, result)
				{
					DataVector myResult(result.getSize());

					/// @todo (heinecke) discuss methods in order to avoid this cast
					((UpDownTwoOpDims*)(this->OpGammaBound))->multParallelBuildingBlock(alpha, myResult, algoDims[i], algoDims[j]);

					// semaphore
					omp_set_lock(&GammaMutex);
					GammaResult.add(myResult);
					omp_unset_lock(&GammaMutex);
				}
			}
		}
	}

	#pragma omp taskwait

	omp_destroy_lock(&GammaMutex);
	omp_destroy_lock(&DeltaMutex);

	// sum up
	result.add(DeltaResult);
	result.sub(GammaResult);
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

void BlackScholesODESolverSystemParallelOMP::mult(DataVector& alpha, DataVector& result)
{
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixInner(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		result.setAll(0.0);

		DataVector temp(result.getSize());
		DataVector temp2(result.getSize());

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
				#pragma omp task shared (alpha, temp)
				{
					applyMassMatrixInner(alpha, temp);
				}

				#pragma omp task shared (alpha, temp2)
				{
					applyLOperatorInner(alpha, temp2);
				}

				#pragma omp taskwait
			}
		}
#else
		applyMassMatrixInner(alpha, temp);
		applyLOperatorInner(alpha, temp2);
#endif

		result.add(temp);
		result.axpy((-1.0)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "CrNic")
	{
		result.setAll(0.0);

		DataVector temp(result.getSize());
		DataVector temp2(result.getSize());

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
				#pragma omp task shared (alpha, temp)
				{
					applyMassMatrixInner(alpha, temp);
				}

				#pragma omp task shared (alpha, temp2)
				{
					applyLOperatorInner(alpha, temp2);
				}

				#pragma omp taskwait
			}
		}
#else
		applyMassMatrixInner(alpha, temp);
		applyLOperatorInner(alpha, temp2);
#endif

		result.add(temp);
		result.axpy((-0.5)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "AdBas")
	{
		result.setAll(0.0);

		applyMassMatrixInner(alpha, result);
	}
	else
	{
		throw new algorithm_exception(" HeatEquationTimestepMatrix::mult : An unknown operation mode was specified!");
	}
}

DataVector* BlackScholesODESolverSystemParallelOMP::generateRHS()
{
	DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(rhs_complete.getSize());
		DataVector temp2(rhs_complete.getSize());

#ifdef USEOMPTHREE
		#pragma omp parallel shared(rhs_complete, temp, temp2)
		{
		#ifndef AIX_XLC
			#pragma omp single nowait
		#endif
		#ifdef AIX_XLC
			#pragma omp single
		#endif
			{
				#pragma omp task shared (temp)
				{
					applyMassMatrixComplete(*this->alpha_complete, temp);
				}

				#pragma omp task shared (temp2)
				{
					applyLOperatorComplete(*this->alpha_complete, temp2);
				}

				#pragma omp taskwait
			}
		}
#else
		applyMassMatrixComplete(*this->alpha_complete, temp);
		applyLOperatorComplete(*this->alpha_complete, temp2);
#endif

		rhs_complete.add(temp);
		rhs_complete.axpy(this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs_complete.setAll(0.0);

		applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(rhs_complete.getSize());
		DataVector temp2(rhs_complete.getSize());

#ifdef USEOMPTHREE
		#pragma omp parallel shared(rhs_complete, temp, temp2)
		{
		#ifndef AIX_XLC
			#pragma omp single nowait
		#endif
		#ifdef AIX_XLC
			#pragma omp single
		#endif
			{
				#pragma omp task shared (temp)
				{
					applyMassMatrixComplete(*this->alpha_complete, temp);
				}

				#pragma omp task shared (temp2)
				{
					applyLOperatorComplete(*this->alpha_complete, temp2);
				}

				#pragma omp taskwait
			}
		}
#else
		applyMassMatrixComplete(*this->alpha_complete, temp);
		applyLOperatorComplete(*this->alpha_complete, temp2);
#endif

		rhs_complete.add(temp);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "AdBas")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);

#ifdef USEOMPTHREE
		#pragma omp parallel shared(rhs_complete, temp)
		{
		#ifndef AIX_XLC
			#pragma omp single nowait
		#endif
		#ifdef AIX_XLC
			#pragma omp single
		#endif
			{
				#pragma omp task shared (temp)
				{
					applyLOperatorComplete(*this->alpha_complete, temp);
				}

				#pragma omp taskwait
			}
		}
#else
		applyLOperatorComplete(*this->alpha_complete, temp);
#endif

		rhs_complete.add(temp);
		temp.mult((2.0)+this->TimestepSize/this->TimestepSize_old);

		DataVector temp_old(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete_old, temp_old);

#ifdef USEOMPTHREE
		#pragma omp parallel shared(rhs_complete, temp_old)
		{
		#ifndef AIX_XLC
			#pragma omp single nowait
		#endif
		#ifdef AIX_XLC
			#pragma omp single
		#endif
			{
				#pragma omp task shared (temp_old)
				{
					applyLOperatorComplete(*this->alpha_complete_old, temp_old);
				}

				#pragma omp taskwait
			}
		}
#else
		applyLOperatorComplete(*this->alpha_complete_old, temp_old);
#endif

		temp_old.mult(this->TimestepSize/this->TimestepSize_old);
		temp.sub(temp_old);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("HeatEquationTimestepMatrix::generateRHS : An unknown operation mode was specified!");
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
		DataVector temp2(alpha_bound.getSize());

#ifdef USEOMPTHREE
		#pragma omp parallel shared(alpha_bound, temp, temp2)
		{
		#ifndef AIX_XLC
			#pragma omp single nowait
		#endif
		#ifdef AIX_XLC
			#pragma omp single
		#endif
			{
				#pragma omp task shared (alpha_bound, temp)
				{
					applyMassMatrixComplete(alpha_bound, temp);
				}

				#pragma omp task shared (alpha_bound, temp2)
				{
					applyLOperatorComplete(alpha_bound, temp2);
				}

				#pragma omp taskwait
			}
		}
#else
		applyMassMatrixComplete(alpha_bound, temp);
		applyLOperatorComplete(alpha_bound, temp2);
#endif
		result_complete.add(temp);
		result_complete.axpy((-1.0)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "CrNic")
	{
		DataVector temp(alpha_bound.getSize());
		DataVector temp2(alpha_bound.getSize());

#ifdef USEOMPTHREE
		#pragma omp parallel shared(alpha_bound, temp, temp2)
		{
		#ifndef AIX_XLC
			#pragma omp single nowait
		#endif
		#ifdef AIX_XLC
			#pragma omp single
		#endif
			{
				#pragma omp task shared (alpha_bound, temp)
				{
					applyMassMatrixComplete(alpha_bound, temp);
				}

				#pragma omp task shared (alpha_bound, temp2)
				{
					applyLOperatorComplete(alpha_bound, temp2);
				}

				#pragma omp taskwait
			}
		}
#else
		applyMassMatrixComplete(alpha_bound, temp);
		applyLOperatorComplete(alpha_bound, temp2);
#endif
		result_complete.add(temp);
		result_complete.axpy((-0.5)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "AdBas")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else
	{
		throw new algorithm_exception("BlackScholesODESolverSystemParallelOMP::generateRHS : An unknown operation mode was specified!");
	}

	rhs_complete.sub(result_complete);

	if (this->rhs != NULL)
	{
		delete this->rhs;
	}

	this->rhs = new DataVector(this->alpha_inner->getSize());
	this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

	return this->rhs;
}

}
