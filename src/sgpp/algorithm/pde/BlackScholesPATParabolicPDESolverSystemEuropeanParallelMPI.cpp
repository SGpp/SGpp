/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/MPI/SGppMPITools.hpp"

#include "algorithm/pde/BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "algorithm/pde/StdUpDown.hpp"
#include "algorithm/pde/UpDownOneOpDim.hpp"
#include "algorithm/pde/UpDownTwoOpDims.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace sg
{
namespace parallel
{

BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& lambda,
			double TimestepSize, std::string OperationMode,
			bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel) : BlackScholesPATParabolicPDESolverSystemEuropean(SparseGrid, alpha, lambda,
			TimestepSize, OperationMode, useCoarsen, coarsenThreshold, adaptSolveMode, numCoarsenPoints, refineThreshold, refineMode, refineMaxLevel)
{}

BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::~BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI()
{
}

void BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	result.setAll(0.0);

	std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();
#ifdef _OPENMP
	omp_lock_t LaplaceMutex;
	omp_init_lock(&LaplaceMutex);
#endif
	sg::base::DataVector LaplaceResult(result);

	// Apply the delta method
	for (size_t i = 0; i < nDims; i++)
	{
		if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
		{
			#pragma omp task firstprivate(i) shared(alpha, LaplaceMutex, LaplaceResult, result, algoDims)
			{
				sg::base::DataVector myResult(result.getSize());

				/// @todo (heinecke) discuss methods in order to avoid this cast
				((sg::pde::UpDownOneOpDim*)(this->OpLaplaceInner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

				// semaphore
#ifdef _OPENMP
				omp_set_lock(&LaplaceMutex);
#endif
				LaplaceResult.add(myResult);
#ifdef _OPENMP
				omp_unset_lock(&LaplaceMutex);
#endif
			}
		}
	}

	#pragma omp taskwait

	result.axpy(-0.5, LaplaceResult);

#ifdef _OPENMP
	omp_destroy_lock(&LaplaceMutex);
#endif
}

void BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	result.setAll(0.0);

	std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
	size_t nDims = algoDims.size();

#ifdef _OPENMP
	omp_lock_t LaplaceMutex;
	omp_init_lock(&LaplaceMutex);
#endif
	sg::base::DataVector LaplaceResult(result);

	// Apply the delta method
	for (size_t i = 0; i < nDims; i++)
	{
		if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
		{
			#pragma omp task firstprivate(i) shared(alpha, LaplaceMutex, LaplaceResult, result, algoDims)
			{
				sg::base::DataVector myResult(result.getSize());

				/// @todo (heinecke) discuss methods in order to avoid this cast
				((sg::pde::UpDownOneOpDim*)(this->OpLaplaceBound))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

				// semaphore
#ifdef _OPENMP
				omp_set_lock(&LaplaceMutex);
#endif
				LaplaceResult.add(myResult);
#ifdef _OPENMP
				omp_unset_lock(&LaplaceMutex);
#endif
			}
		}
	}

	#pragma omp taskwait

	result.axpy(-0.5, LaplaceResult);

#ifdef _OPENMP
	omp_destroy_lock(&LaplaceMutex);
#endif
}

void BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());
	result.setAll(0.0);
	size_t nDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions().size();

	if (nDims % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
	{
		((sg::pde::StdUpDown*)(this->OpLTwoInner))->multParallelBuildingBlock(alpha, temp);

		result.add(temp);
	}
}

void BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());
	result.setAll(0.0);
	size_t nDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions().size();

	if (nDims % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank())
	{
		((sg::pde::StdUpDown*)(this->OpLTwoBound))->multParallelBuildingBlock(alpha, temp);

		result.add(temp);
	}
}

void BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	// distribute the current grid coefficients
	myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);

	result.setAll(0.0);

	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixInner(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		sg::base::DataVector temp(result.getSize());
		sg::base::DataVector temp2(result.getSize());

		#pragma omp parallel shared(alpha, result)
		{
			#pragma omp single nowait
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

		result.add(temp);
		result.axpy((-1.0)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "CrNic")
	{
		sg::base::DataVector temp(result.getSize());
		sg::base::DataVector temp2(result.getSize());

		#pragma omp parallel shared(alpha, result)
		{
			#pragma omp single nowait
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
		throw new sg::base::algorithm_exception(" BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::mult : An unknown operation mode was specified!");
	}

	// aggregate all results
	myGlobalMPIComm->reduceGridCoefficientsOnRank0(result);
}

sg::base::DataVector* BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::generateRHS()
{
	// distribute the current grid coefficients
	myGlobalMPIComm->broadcastGridCoefficientsFromRank0(*(this->alpha_complete));

	sg::base::DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		sg::base::DataVector temp(rhs_complete.getSize());
		sg::base::DataVector temp2(rhs_complete.getSize());
		sg::base::DataVector myAlpha(*this->alpha_complete);

		#pragma omp parallel shared(myAlpha, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (myAlpha, temp)
				{
					applyMassMatrixComplete(myAlpha, temp);
				}

				#pragma omp task shared (myAlpha, temp2)
				{
					applyLOperatorComplete(myAlpha, temp2);
				}

				#pragma omp taskwait
			}
		}

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

		sg::base::DataVector temp(rhs_complete.getSize());
		sg::base::DataVector temp2(rhs_complete.getSize());
		sg::base::DataVector myAlpha(*this->alpha_complete);

		#pragma omp parallel shared(myAlpha, temp, temp2)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (myAlpha, temp)
				{
					applyMassMatrixComplete(myAlpha, temp);
				}

				#pragma omp task shared (myAlpha, temp2)
				{
					applyLOperatorComplete(myAlpha, temp2);
				}

				#pragma omp taskwait
			}
		}

		rhs_complete.add(temp);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "AdBas")
	{
		rhs_complete.setAll(0.0);

		sg::base::DataVector temp(this->alpha_complete->getSize());
		sg::base::DataVector myAlpha(*this->alpha_complete);
		sg::base::DataVector myOldAlpha(*this->alpha_complete_old);

		applyMassMatrixComplete(*this->alpha_complete, temp);

		#pragma omp parallel shared(myAlpha, temp)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (myAlpha, temp)
				{
					applyLOperatorComplete(myAlpha, temp);
				}

				#pragma omp taskwait
			}
		}

		rhs_complete.add(temp);
		temp.mult((2.0)+this->TimestepSize/this->TimestepSize_old);

		sg::base::DataVector temp_old(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete_old, temp_old);

		#pragma omp parallel shared(myOldAlpha, temp_old)
		{
			#pragma omp single nowait
			{
				#pragma omp task shared (myOldAlpha, temp_old)
				{
					applyLOperatorComplete(myOldAlpha, temp_old);
				}

				#pragma omp taskwait
			}
		}

		temp_old.mult(this->TimestepSize/this->TimestepSize_old);
		temp.sub(temp_old);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new sg::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::generateRHS : An unknown operation mode was specified!");
	}

	// aggregate all results
	myGlobalMPIComm->reduceGridCoefficientsOnRank0(rhs_complete);

	// Now we have the right hand side, lets apply the riskfree rate for the next timestep
	this->startTimestep();

	// Now apply the boundary ansatzfunctions to the inner ansatzfunctions
	sg::base::DataVector result_complete(this->alpha_complete->getSize());
	sg::base::DataVector alpha_bound(*this->alpha_complete);

	result_complete.setAll(0.0);

	this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

	// apply CG Matrix
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else if (this->tOperationMode == "ImEul")
	{
		sg::base::DataVector temp(alpha_bound.getSize());
		sg::base::DataVector temp2(alpha_bound.getSize());

		#pragma omp parallel shared(alpha_bound, temp, temp2)
		{
			#pragma omp single nowait
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

		result_complete.add(temp);
		result_complete.axpy((-1.0)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "CrNic")
	{
		sg::base::DataVector temp(alpha_bound.getSize());
		sg::base::DataVector temp2(alpha_bound.getSize());

		#pragma omp parallel shared(alpha_bound, temp, temp2)
		{
			#pragma omp single nowait
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

		result_complete.add(temp);
		result_complete.axpy((-0.5)*this->TimestepSize, temp2);
	}
	else if (this->tOperationMode == "AdBas")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else
	{
		throw new sg::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::generateRHS : An unknown operation mode was specified!");
	}

	// aggregate all results
	myGlobalMPIComm->reduceGridCoefficientsOnRank0(result_complete);

	rhs_complete.sub(result_complete);

	if (this->rhs != NULL)
	{
		delete this->rhs;
	}

	this->rhs = new sg::base::DataVector(this->alpha_inner->getSize());

	if (myGlobalMPIComm->getMyRank() == 0)
	{
		this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);
	}
	else
	{
		this->rhs->setAll(0.0);
	}

	return this->rhs;
}

void BlackScholesPATParabolicPDESolverSystemEuropeanParallelMPI::finishTimestep(bool isLastTimestep)
{
	// Adaptivity stuff is done on rank 0 only
	if (myGlobalMPIComm->getMyRank() == 0)
	{
		// Replace the inner coefficients on the boundary grid
		this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

		// add number of Gridpoints
		this->numSumGridpointsInner += this->InnerGrid->getSize();
		this->numSumGridpointsComplete += this->BoundGrid->getSize();

		if (this->useCoarsen == true && isLastTimestep == false)
		{
			///////////////////////////////////////////////////
			// Start integrated refinement & coarsening
			///////////////////////////////////////////////////

			size_t originalGridSize = this->BoundGrid->getStorage()->size();

			// Coarsen the grid
			base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

			//std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
			//std::cout << "Grid Size: " << originalGridSize << std::endl;

			if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine")
			{
				size_t numRefines = myGenerator->getNumberOfRefinablePoints();
				base::SurplusRefinementFunctor* myRefineFunc = new base::SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refineThreshold);
				if (this->refineMode == "maxLevel")
				{
					myGenerator->refineMaxLevel(myRefineFunc, this->refineMaxLevel);
					this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
				}
				if (this->refineMode == "classic")
				{
					myGenerator->refine(myRefineFunc);
					this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
				}
				delete myRefineFunc;
			}

			if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine")
			{
				size_t numCoarsen = myGenerator->getNumberOfRemoveablePoints();
				base::SurplusCoarseningFunctor* myCoarsenFunctor = new base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);
				myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
				delete myCoarsenFunctor;
			}

			delete myGenerator;

			///////////////////////////////////////////////////
			// End integrated refinement & coarsening
			///////////////////////////////////////////////////
		}
	}

	// only communicate if adaptivity is switched on
	if (this->useCoarsen == true && isLastTimestep == false)
	{
		// Communicate new grid
		if (myGlobalMPIComm->getMyRank() == 0)
		{
			std::string bound_grid_storage = this->BoundGrid->getStorage()->serialize();

			myGlobalMPIComm->broadcastGridStorage(bound_grid_storage);

			// rebuild the inner grid + coefficients
			this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
		}
		else
		{
			std::string bound_grid_storage = "";

			myGlobalMPIComm->receiveGridStorage(bound_grid_storage);

			this->BoundGrid->getStorage()->emptyStorage();
			this->BoundGrid->getStorage()->unserialize_noAlgoDims(bound_grid_storage);
			this->alpha_complete->resize(this->BoundGrid->getStorage()->size());

			// rebuild the inner grid + coefficients
			this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
		}

		myGlobalMPIComm->Barrier();
	}
}

}

}
