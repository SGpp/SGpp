// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/pde/algorithm/StdUpDown.hpp>
#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>
#ifdef USE_ENHANCED_UPDOWN
#include <sgpp/misc/pde/algorithm/UpDownOneOpDimEnhanced.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>

namespace SGPP {
namespace finance {

BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::
    BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP(
        SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& lambda,
        SGPP::base::DataMatrix& eigenvecs, SGPP::base::DataVector& mu_hat, float_t TimestepSize,
        std::string OperationMode, float_t dStrike, std::string option_type, float_t r,
        bool useCoarsen, float_t coarsenThreshold, std::string adaptSolveMode, int numCoarsenPoints,
        float_t refineThreshold, std::string refineMode,
        SGPP::base::GridIndex::level_type refineMaxLevel)
    : BlackScholesPATParabolicPDESolverSystemEuroAmer(
          SparseGrid, alpha, lambda, eigenvecs, mu_hat, TimestepSize, OperationMode, dStrike,
          option_type, r, useCoarsen, coarsenThreshold, adaptSolveMode, numCoarsenPoints,
          refineThreshold, refineMode, refineMaxLevel),
      rhs_corrector(NULL) {}

BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::
    ~BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP() {
  if (this->rhs_corrector != NULL) {
    delete this->rhs_corrector;
  }
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::applyLOperatorInner(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  result.setAll(0.0);

#ifdef USE_ENHANCED_UPDOWN
  ((SGPP::pde::UpDownOneOpDimEnhanced*)(this->OpLaplaceInner))
      ->multParallelBuildingBlock(alpha, result);
  result.mult(-0.5);
#else
  std::vector<size_t> algoDims = this->InnerGrid->getStorage().getAlgorithmicDimensions();
  size_t nDims = algoDims.size();

#ifdef _OPENMP
  omp_lock_t LaplaceMutex;
  omp_init_lock(&LaplaceMutex);
#endif
  SGPP::base::DataVector LaplaceResult(result);

  // Apply Laplace
  for (size_t i = 0; i < nDims; i++) {
#pragma omp task firstprivate(i) shared(alpha, LaplaceMutex, LaplaceResult, result, algoDims)
    {
      SGPP::base::DataVector myResult(result.getSize());

      /// discuss methods in order to avoid this cast
      ((SGPP::pde::UpDownOneOpDim*)(this->OpLaplaceInner))
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

#pragma omp taskwait

  result.axpy(-0.5, LaplaceResult);
#ifdef _OPENMP
  omp_destroy_lock(&LaplaceMutex);
#endif
#endif
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::applyLOperatorComplete(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  result.setAll(0.0);

#ifdef USE_ENHANCED_UPDOWN
  ((SGPP::pde::UpDownOneOpDimEnhanced*)(this->OpLaplaceBound))
      ->multParallelBuildingBlock(alpha, result);
  result.mult(-0.5);
#else
  std::vector<size_t> algoDims = this->BoundGrid->getStorage().getAlgorithmicDimensions();
  size_t nDims = algoDims.size();
#ifdef _OPENMP
  omp_lock_t LaplaceMutex;
  omp_init_lock(&LaplaceMutex);
#endif
  SGPP::base::DataVector LaplaceResult(result);

  // Apply Laplace
  for (size_t i = 0; i < nDims; i++) {
#pragma omp task firstprivate(i) shared(alpha, LaplaceMutex, LaplaceResult, result, algoDims)
    {
      SGPP::base::DataVector myResult(result.getSize());

      /// discuss methods in order to avoid this cast
      ((SGPP::pde::UpDownOneOpDim*)(this->OpLaplaceBound))
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

#pragma omp taskwait

  result.axpy(-0.5, LaplaceResult);

#ifdef _OPENMP
  omp_destroy_lock(&LaplaceMutex);
#endif
#endif
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::applyMassMatrixInner(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  ((SGPP::pde::StdUpDown*)(this->OpLTwoInner))->multParallelBuildingBlock(alpha, temp);

  result.add(temp);
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::applyMassMatrixComplete(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  ((SGPP::pde::StdUpDown*)(this->OpLTwoBound))->multParallelBuildingBlock(alpha, temp);

  result.add(temp);
}

void BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::mult(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  result.setAll(0.0);

  if (this->tOperationMode == "ExEul") {
    applyMassMatrixInner(alpha, result);
  } else if (this->tOperationMode == "ImEul") {
    SGPP::base::DataVector temp(result.getSize());
    SGPP::base::DataVector temp2(result.getSize());

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
    SGPP::base::DataVector temp(result.getSize());
    SGPP::base::DataVector temp2(result.getSize());

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
    throw new SGPP::base::algorithm_exception(
        " BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::mult : An unknown operation "
        "mode was specified!");
  }
}

SGPP::base::DataVector* BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP::generateRHS() {
  SGPP::base::DataVector rhs_complete(this->alpha_complete->getSize());

  if (this->tOperationMode == "ExEul") {
    rhs_complete.setAll(0.0);

    SGPP::base::DataVector temp(rhs_complete.getSize());
    SGPP::base::DataVector temp2(rhs_complete.getSize());
    SGPP::base::DataVector myAlpha(*this->alpha_complete);

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

    SGPP::base::DataVector myAlpha(*this->alpha_complete);

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

    SGPP::base::DataVector temp(rhs_complete.getSize());
    SGPP::base::DataVector temp2(rhs_complete.getSize());
    SGPP::base::DataVector myAlpha(*this->alpha_complete);

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

    SGPP::base::DataVector temp(this->alpha_complete->getSize());
    SGPP::base::DataVector myAlpha(*this->alpha_complete);
    SGPP::base::DataVector myOldAlpha(*this->alpha_complete_old);

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

    SGPP::base::DataVector temp_old(this->alpha_complete->getSize());

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
    throw new SGPP::base::algorithm_exception(
        "BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::generateRHS : An unknown "
        "operation mode was specified!");
  }

  if (this->useCoarsen == true || this->nExecTimesteps == 0 || this->bnewODESolver == true ||
      this->option_type == "std_amer_put") {
    this->bnewODESolver = false;

    // Now apply the boundary ansatzfunctions to the inner ansatzfunctions
    SGPP::base::DataVector result_complete(this->alpha_complete->getSize());
    SGPP::base::DataVector alpha_bound(*this->alpha_complete);

    result_complete.setAll(0.0);

    this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

    // apply CG Matrix
    if (this->tOperationMode == "ExEul") {
      applyMassMatrixComplete(alpha_bound, result_complete);
    } else if (this->tOperationMode == "ImEul") {
      SGPP::base::DataVector temp(alpha_bound.getSize());
      SGPP::base::DataVector temp2(alpha_bound.getSize());

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
      SGPP::base::DataVector temp(alpha_bound.getSize());
      SGPP::base::DataVector temp2(alpha_bound.getSize());

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
      throw new SGPP::base::algorithm_exception(
          "BlackScholesPATParabolicPDESolverSystemEuropeanParallelOMP::generateRHS : An unknown "
          "operation mode was specified!");
    }

    // Store right hand side corrector
    if (this->rhs_corrector != NULL) {
      delete this->rhs_corrector;
    }

    this->rhs_corrector = new SGPP::base::DataVector(this->alpha_complete->getSize());
    *(this->rhs_corrector) = result_complete;
  }

  rhs_complete.sub(*(this->rhs_corrector));

  //  rhs_complete.sub(result_complete);

  if (this->rhs != NULL) {
    delete this->rhs;
  }

  this->rhs = new SGPP::base::DataVector(this->alpha_inner->getSize());
  this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

  return this->rhs;
}
}  // namespace finance
}  // namespace SGPP
