// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>

#include <sgpp/parallel/finance/algorithm/BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/pde/algorithm/StdUpDown.hpp>
#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>
#include <sgpp/pde/algorithm/UpDownTwoOpDims.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <algorithm>
#include <vector>

namespace sgpp {
namespace parallel {

BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::
    BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI(
        sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& alpha, sgpp::base::DataVector& lambda,
        sgpp::base::DataMatrix& eigenvecs, sgpp::base::DataVector& mu_hat, double TimestepSize,
        std::string OperationMode, double dStrike, std::string option_type, double r,
        bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode, int numCoarsenPoints,
        double refineThreshold, std::string refineMode,
        sgpp::base::GridPoint::level_type refineMaxLevel)
    : BlackScholesPATParabolicPDESolverSystemEuroAmer(
          SparseGrid, alpha, lambda, eigenvecs, mu_hat, TimestepSize, OperationMode, dStrike,
          option_type, r, useCoarsen, coarsenThreshold, adaptSolveMode, numCoarsenPoints,
          refineThreshold, refineMode, refineMaxLevel),
      rhs_corrector(NULL) {}

BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::
    ~BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI() {
  if (this->rhs_corrector != NULL) {
    delete this->rhs_corrector;
  }
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::applyLOperatorInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

  std::vector<size_t> algoDims = this->InnerGrid->getStorage().getAlgorithmicDimensions();
  size_t nDims = algoDims.size();
#ifdef _OPENMP
  omp_lock_t LaplaceMutex;
  omp_init_lock(&LaplaceMutex);
#endif
  sgpp::base::DataVector LaplaceResult(result);

  // Apply the delta method
  for (size_t i = 0; i < nDims; i++) {
    if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
#pragma omp task firstprivate(i) shared(alpha, LaplaceMutex, LaplaceResult, result, algoDims)
      {
        sgpp::base::DataVector myResult(result.getSize());

        /// discuss methods in order to avoid this cast
        ((sgpp::pde::UpDownOneOpDim*)(this->OpLaplaceInner))
            ->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

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

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::applyLOperatorComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

  std::vector<size_t> algoDims = this->InnerGrid->getStorage().getAlgorithmicDimensions();
  size_t nDims = algoDims.size();

#ifdef _OPENMP
  omp_lock_t LaplaceMutex;
  omp_init_lock(&LaplaceMutex);
#endif
  sgpp::base::DataVector LaplaceResult(result);

  // Apply the delta method
  for (size_t i = 0; i < nDims; i++) {
    if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
#pragma omp task firstprivate(i) shared(alpha, LaplaceMutex, LaplaceResult, result, algoDims)
      {
        sgpp::base::DataVector myResult(result.getSize());

        /// discuss methods in order to avoid this cast
        ((sgpp::pde::UpDownOneOpDim*)(this->OpLaplaceBound))
            ->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

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

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::applyMassMatrixInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());
  result.setAll(0.0);
  size_t nDims = this->InnerGrid->getStorage().getAlgorithmicDimensions().size();

  if (nDims % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
    ((sgpp::pde::StdUpDown*)(this->OpLTwoInner))->multParallelBuildingBlock(alpha, temp);

    result.add(temp);
  }
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::applyMassMatrixComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());
  result.setAll(0.0);
  size_t nDims = this->InnerGrid->getStorage().getAlgorithmicDimensions().size();

  if (nDims % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
    ((sgpp::pde::StdUpDown*)(this->OpLTwoBound))->multParallelBuildingBlock(alpha, temp);

    result.add(temp);
  }
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::mult(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  // distribute the current grid coefficients
  // myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);

  result.setAll(0.0);

  if (this->tOperationMode == "ExEul") {
    applyMassMatrixInner(alpha, result);
  } else if (this->tOperationMode == "ImEul") {
    sgpp::base::DataVector temp(result.getSize());
    sgpp::base::DataVector temp2(result.getSize());

#pragma omp parallel shared(alpha, result)
    {
#pragma omp single nowait
      {
#pragma omp task shared(alpha, temp)
        { applyMassMatrixInner(alpha, temp); }

#pragma omp task shared(alpha, temp2)
        { applyLOperatorInner(alpha, temp2); }

#pragma omp taskwait
      }
    }

    result.add(temp);
    result.axpy((-1.0) * this->TimestepSize, temp2);
  } else if (this->tOperationMode == "CrNic") {
    sgpp::base::DataVector temp(result.getSize());
    sgpp::base::DataVector temp2(result.getSize());

#pragma omp parallel shared(alpha, result)
    {
#pragma omp single nowait
      {
#pragma omp task shared(alpha, temp)
        { applyMassMatrixInner(alpha, temp); }

#pragma omp task shared(alpha, temp2)
        { applyLOperatorInner(alpha, temp2); }

#pragma omp taskwait
      }
    }

    result.add(temp);
    result.axpy((-0.5) * this->TimestepSize, temp2);
  } else if (this->tOperationMode == "AdBas") {
    result.setAll(0.0);

    applyMassMatrixInner(alpha, result);
  } else {
    throw sgpp::base::algorithm_exception(
        " BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::mult : An unknown operation "
        "mode was specified!");
  }

  // aggregate all results
  // myGlobalMPIComm->reduceGridCoefficientsOnRank0(result);
  myGlobalMPIComm->reduceGridCoefficients(result);
}

sgpp::base::DataVector* BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::generateRHS() {
  // We have no adaptivity during solving and we already have a rhs
  // --> just return
  // if (this->rhs != NULL && this->useCoarsen == false)
  //{
  //  return this->rhs;
  //}

  // distribute the current grid coefficients
  myGlobalMPIComm->broadcastGridCoefficientsFromRank0(*(this->alpha_complete));

  sgpp::base::DataVector rhs_complete(this->alpha_complete->getSize());

  if (this->tOperationMode == "ExEul") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(rhs_complete.getSize());
    sgpp::base::DataVector temp2(rhs_complete.getSize());
    sgpp::base::DataVector myAlpha(*this->alpha_complete);

#pragma omp parallel shared(myAlpha, temp, temp2)
    {
#pragma omp single nowait
      {
#pragma omp task shared(myAlpha, temp)
        { applyMassMatrixComplete(myAlpha, temp); }

#pragma omp task shared(myAlpha, temp2)
        { applyLOperatorComplete(myAlpha, temp2); }

#pragma omp taskwait
      }
    }

    rhs_complete.add(temp);
    rhs_complete.axpy(this->TimestepSize, temp2);
  } else if (this->tOperationMode == "ImEul") {
    rhs_complete.setAll(0.0);
    sgpp::base::DataVector myAlpha(*this->alpha_complete);

#pragma omp parallel shared(myAlpha, rhs_complete)
    {
#pragma omp single nowait
      {
#pragma omp task shared(rhs_complete)
        { applyMassMatrixComplete(myAlpha, rhs_complete); }

#pragma omp taskwait
      }
    }
  } else if (this->tOperationMode == "CrNic") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(rhs_complete.getSize());
    sgpp::base::DataVector temp2(rhs_complete.getSize());
    sgpp::base::DataVector myAlpha(*this->alpha_complete);

#pragma omp parallel shared(myAlpha, temp, temp2)
    {
#pragma omp single nowait
      {
#pragma omp task shared(myAlpha, temp)
        { applyMassMatrixComplete(myAlpha, temp); }

#pragma omp task shared(myAlpha, temp2)
        { applyLOperatorComplete(myAlpha, temp2); }

#pragma omp taskwait
      }
    }

    rhs_complete.add(temp);
    rhs_complete.axpy((0.5) * this->TimestepSize, temp2);
  } else if (this->tOperationMode == "AdBas") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());
    sgpp::base::DataVector myAlpha(*this->alpha_complete);
    sgpp::base::DataVector myOldAlpha(*this->alpha_complete_old);

    applyMassMatrixComplete(*this->alpha_complete, temp);

#pragma omp parallel shared(myAlpha, temp)
    {
#pragma omp single nowait
      {
#pragma omp task shared(myAlpha, temp)
        { applyLOperatorComplete(myAlpha, temp); }

#pragma omp taskwait
      }
    }

    rhs_complete.add(temp);
    temp.mult((2.0) + this->TimestepSize / this->TimestepSize_old);

    sgpp::base::DataVector temp_old(this->alpha_complete->getSize());

    applyMassMatrixComplete(*this->alpha_complete_old, temp_old);

#pragma omp parallel shared(myOldAlpha, temp_old)
    {
#pragma omp single nowait
      {
#pragma omp task shared(myOldAlpha, temp_old)
        { applyLOperatorComplete(myOldAlpha, temp_old); }

#pragma omp taskwait
      }
    }

    temp_old.mult(this->TimestepSize / this->TimestepSize_old);
    temp.sub(temp_old);
    rhs_complete.axpy((0.5) * this->TimestepSize, temp);
  } else {
    throw sgpp::base::algorithm_exception(
        "BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::generateRHS : An unknown "
        "operation mode was specified!");
  }

  // aggregate all results
  myGlobalMPIComm->reduceGridCoefficientsOnRank0(rhs_complete);

  if (this->useCoarsen == true || this->nExecTimesteps == 0 || this->bnewODESolver == true ||
      this->option_type == "std_amer_put") {
    this->bnewODESolver = false;

    // Now apply the boundary ansatzfunctions to the inner ansatzfunctions
    sgpp::base::DataVector result_complete(this->alpha_complete->getSize());
    sgpp::base::DataVector alpha_bound(*this->alpha_complete);

    result_complete.setAll(0.0);

    this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

    // apply CG Matrix
    if (this->tOperationMode == "ExEul") {
      applyMassMatrixComplete(alpha_bound, result_complete);
    } else if (this->tOperationMode == "ImEul") {
      sgpp::base::DataVector temp(alpha_bound.getSize());
      sgpp::base::DataVector temp2(alpha_bound.getSize());

#pragma omp parallel shared(alpha_bound, temp, temp2)
      {
#pragma omp single nowait
        {
#pragma omp task shared(alpha_bound, temp)
          { applyMassMatrixComplete(alpha_bound, temp); }

#pragma omp task shared(alpha_bound, temp2)
          { applyLOperatorComplete(alpha_bound, temp2); }

#pragma omp taskwait
        }
      }

      result_complete.add(temp);
      result_complete.axpy((-1.0) * this->TimestepSize, temp2);
    } else if (this->tOperationMode == "CrNic") {
      sgpp::base::DataVector temp(alpha_bound.getSize());
      sgpp::base::DataVector temp2(alpha_bound.getSize());

#pragma omp parallel shared(alpha_bound, temp, temp2)
      {
#pragma omp single nowait
        {
#pragma omp task shared(alpha_bound, temp)
          { applyMassMatrixComplete(alpha_bound, temp); }

#pragma omp task shared(alpha_bound, temp2)
          { applyLOperatorComplete(alpha_bound, temp2); }

#pragma omp taskwait
        }
      }

      result_complete.add(temp);
      result_complete.axpy((-0.5) * this->TimestepSize, temp2);
    } else if (this->tOperationMode == "AdBas") {
      applyMassMatrixComplete(alpha_bound, result_complete);
    } else {
      throw sgpp::base::algorithm_exception(
          "BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::generateRHS : An unknown "
          "operation mode was specified!");
    }

    // aggregate all results
    myGlobalMPIComm->reduceGridCoefficientsOnRank0(result_complete);

    // Store right hand side corrector
    if (this->rhs_corrector != NULL) {
      delete this->rhs_corrector;
    }

    this->rhs_corrector = new sgpp::base::DataVector(this->alpha_complete->getSize());
    *(this->rhs_corrector) = result_complete;
  }

  rhs_complete.sub(*(this->rhs_corrector));

  //  rhs_complete.sub(result_complete);

  if (this->rhs != NULL) {
    delete this->rhs;
  }

  this->rhs = new sgpp::base::DataVector(this->alpha_inner->getSize());

  if (myGlobalMPIComm->getMyRank() == 0) {
    this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);
  } else {
    this->rhs->setAll(0.0);
  }

  return this->rhs;
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI::finishTimestep(
    bool isLastTimestep) {
  this->nExecTimesteps++;

  // Adaptivity stuff is done on rank 0 only
  if (myGlobalMPIComm->getMyRank() == 0) {
    // Replace the inner coefficients on the boundary grid
    this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

    // check if we are doing an American put -> handle early exercise
    if (this->option_type == "std_amer_put") {
      double current_time = static_cast<double>(this->nExecTimesteps) * this->TimestepSize;

      std::unique_ptr<sgpp::base::OperationHierarchisation> myHierarchisation(
          sgpp::op_factory::createOperationHierarchisation(*this->BoundGrid));
      myHierarchisation->doDehierarchisation(*this->alpha_complete);
      size_t dim = this->BoundGrid->getStorage().getDimension();

      double* coords_val = new double[dim];

      for (size_t i = 0; i < this->BoundGrid->getStorage().getSize(); i++) {
        std::vector<double> eval_point_coord;
        std::string coords = this->BoundGrid->getStorage().getCoordinates(
            this->BoundGrid->getStorage().getPoint(i)).toString();
        std::stringstream coordsStream(coords);

        double tmp;

        // read coordinates
        for (size_t j = 0; j < dim; j++) {
          coordsStream >> tmp;

          coords_val[j] = tmp;
        }

        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          double inner_tmp = 0.0;

          for (size_t l = 0; l < dim; l++) {
            inner_tmp += this->eigenvecs->get(j, l) *
                         (coords_val[l] - (current_time * this->mu_hat->get(l)));
          }

          tmp += exp(inner_tmp);
        }

        double payoff = std::max<double>(this->dStrike - (tmp / static_cast<double>(dim)),
                                         0.0);  // *exp(((-1.0)*(this->r*current_time)));
        double discounted_value =
            ((*this->alpha_complete)[i]) * exp(((-1.0) * (this->r * this->TimestepSize)));

        (*this->alpha_complete)[i] = std::max<double>(payoff, discounted_value);
      }

      delete[] coords_val;

      myHierarchisation->doHierarchisation(*this->alpha_complete);
    }

    // add number of Gridpoints
    this->numSumGridpointsInner += this->InnerGrid->getSize();
    this->numSumGridpointsComplete += this->BoundGrid->getSize();

    if (this->useCoarsen == true && isLastTimestep == false) {
      ///////////////////////////////////////////////////
      // Start integrated refinement & coarsening
      ///////////////////////////////////////////////////

      size_t originalGridSize = this->BoundGrid->getStorage().getSize();

      // Coarsen the grid
      base::GridGenerator& myGenerator = this->BoundGrid->getGenerator();

      // std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
      // std::cout << "Grid Size: " << originalGridSize << std::endl;

      if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine") {
        size_t numRefines = myGenerator.getNumberOfRefinablePoints();
        base::SurplusRefinementFunctor myRefineFunc(*this->alpha_complete, numRefines,
                                                    this->refineThreshold);

        if (this->refineMode == "maxLevel") {
          myGenerator.refineMaxLevel(myRefineFunc, this->refineMaxLevel);
          this->alpha_complete->resizeZero(this->BoundGrid->getStorage().getSize());
        }

        if (this->refineMode == "classic") {
          myGenerator.refine(myRefineFunc);
          this->alpha_complete->resizeZero(this->BoundGrid->getStorage().getSize());
        }
      }

      if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine") {
        size_t numCoarsen = myGenerator.getNumberOfRemovablePoints();
        base::SurplusCoarseningFunctor myCoarsenFunctor(*this->alpha_complete, numCoarsen,
                                                        this->coarsenThreshold);
        myGenerator.coarsenNFirstOnly(myCoarsenFunctor, *this->alpha_complete, originalGridSize);
      }

      ///////////////////////////////////////////////////
      // End integrated refinement & coarsening
      ///////////////////////////////////////////////////
    }
  }

  // only communicate if adaptivity is switched on
  if (this->useCoarsen == true && isLastTimestep == false) {
    // Communicate new grid
    if (myGlobalMPIComm->getMyRank() == 0) {
      std::string bound_grid_storage = this->BoundGrid->getStorage().serialize();

      myGlobalMPIComm->broadcastGridStorage(bound_grid_storage);

      // rebuild the inner grid + coefficients
      this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                                     &this->InnerGrid, &this->alpha_inner);
    } else {
      std::string bound_grid_storage = "";

      myGlobalMPIComm->receiveGridStorage(bound_grid_storage);

      this->BoundGrid->getStorage().clear();
      this->BoundGrid->getStorage().unserializeNoAlgoDims(bound_grid_storage);
      this->alpha_complete->resize(this->BoundGrid->getStorage().getSize());

      // rebuild the inner grid + coefficients
      this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                                     &this->InnerGrid, &this->alpha_inner);
    }
  }
}
}  // namespace parallel
}  // namespace sgpp
