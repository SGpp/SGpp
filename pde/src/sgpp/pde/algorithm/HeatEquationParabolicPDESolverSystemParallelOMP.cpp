// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystemParallelOMP.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/pde/algorithm/StdUpDown.hpp>
#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/pde/operation/PdeOpFactory.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sgpp/globaldef.hpp>

#include <vector>
#include <string>

namespace sgpp {
namespace pde {

HeatEquationParabolicPDESolverSystemParallelOMP::HeatEquationParabolicPDESolverSystemParallelOMP(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& alpha, double a, double TimestepSize,
    std::string OperationMode) {
  this->a = a;
  this->tOperationMode = OperationMode;
  this->TimestepSize = TimestepSize;
  this->BoundGrid = &SparseGrid;
  this->alpha_complete = &alpha;
  this->InnerGrid = nullptr;
  this->alpha_inner = nullptr;

  this->BoundaryUpdate = new sgpp::base::DirichletUpdateVector(SparseGrid.getStorage());
  this->GridConverter = new sgpp::base::DirichletGridConverter();

  this->OpLaplaceBound = op_factory::createOperationLaplace(SparseGrid);
  this->OpMassBound = sgpp::op_factory::createOperationLTwoDotProduct(SparseGrid);

  // create the inner grid
  this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                               &this->InnerGrid, &this->alpha_inner);

  // Create needed operations, on inner grid
  this->OpLaplaceInner = op_factory::createOperationLaplace(*this->InnerGrid);
  this->OpMassInner = sgpp::op_factory::createOperationLTwoDotProduct(*this->InnerGrid);

  // right hand side if System
  this->rhs = nullptr;
}

HeatEquationParabolicPDESolverSystemParallelOMP::
    ~HeatEquationParabolicPDESolverSystemParallelOMP() {
  delete this->OpLaplaceBound;
  delete this->OpMassBound;
  delete this->OpLaplaceInner;
  delete this->OpMassInner;

  delete this->BoundaryUpdate;
  delete this->GridConverter;

  if (this->InnerGrid != nullptr) {
    delete this->InnerGrid;
  }

  if (this->alpha_inner != nullptr) {
    delete this->alpha_inner;
  }

  if (this->rhs != nullptr) {
    delete this->rhs;
  }
}

void HeatEquationParabolicPDESolverSystemParallelOMP::applyMassMatrixComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  reinterpret_cast<StdUpDown*>(this->OpMassBound)->multParallelBuildingBlock(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystemParallelOMP::applyLOperatorComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  std::vector<size_t> algoDims = this->InnerGrid->getStorage().getAlgorithmicDimensions();
  size_t nDims = algoDims.size();
#ifdef _OPENMP
  omp_lock_t Mutex;
  omp_init_lock(&Mutex);
#endif

  // Apply Laplace, parallel in Dimensions
  for (size_t i = 0; i < nDims; i++) {
#pragma omp task firstprivate(i) shared(alpha, temp, result, algoDims)
    {
      sgpp::base::DataVector myResult(result.getSize());

      /// discuss methods in order to avoid this cast
      reinterpret_cast<UpDownOneOpDim*>(this->OpLaplaceBound)
          ->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

// semaphore
#ifdef _OPENMP
      omp_set_lock(&Mutex);
#endif
      temp.add(myResult);
#ifdef _OPENMP
      omp_unset_lock(&Mutex);
#endif
    }
  }

#pragma omp taskwait

#ifdef _OPENMP
  omp_destroy_lock(&Mutex);
#endif

  result.axpy((-1.0) * this->a, temp);
}

void HeatEquationParabolicPDESolverSystemParallelOMP::applyMassMatrixInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  reinterpret_cast<StdUpDown*>(this->OpMassInner)->multParallelBuildingBlock(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystemParallelOMP::applyLOperatorInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  std::vector<size_t> algoDims = this->InnerGrid->getStorage().getAlgorithmicDimensions();
  size_t nDims = algoDims.size();
#ifdef _OPENMP
  omp_lock_t Mutex;
  omp_init_lock(&Mutex);
#endif

  // Apply Laplace, parallel in Dimensions
  for (size_t i = 0; i < nDims; i++) {
#pragma omp task firstprivate(i) shared(alpha, temp, result, algoDims)
    {
      sgpp::base::DataVector myResult(result.getSize());

      /// discuss methods in order to avoid this cast
      reinterpret_cast<UpDownOneOpDim*>(this->OpLaplaceInner)
          ->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

// semaphore
#ifdef _OPENMP
      omp_set_lock(&Mutex);
#endif
      temp.add(myResult);
#ifdef _OPENMP
      omp_unset_lock(&Mutex);
#endif
    }
  }

#pragma omp taskwait

#ifdef _OPENMP
  omp_destroy_lock(&Mutex);
#endif

  result.axpy((-1.0) * this->a, temp);
}

void HeatEquationParabolicPDESolverSystemParallelOMP::finishTimestep() {
  // Replace the inner coefficients on the boundary grid
  this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void HeatEquationParabolicPDESolverSystemParallelOMP::coarsenAndRefine(bool isLastTimestep) {}

void HeatEquationParabolicPDESolverSystemParallelOMP::startTimestep() {}

void HeatEquationParabolicPDESolverSystemParallelOMP::mult(sgpp::base::DataVector& alpha,
                                                           sgpp::base::DataVector& result) {
  result.setAll(0.0);

  if (this->tOperationMode == "ExEul") {
    applyMassMatrixInner(alpha, result);
  } else if (this->tOperationMode == "ImEul") {
    sgpp::base::DataVector temp(result.getSize());
    sgpp::base::DataVector temp2(result.getSize());

#pragma omp parallel shared(alpha, result, temp, temp2)
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

#pragma omp parallel shared(alpha, result, temp, temp2)
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
  } else {
    throw sgpp::base::algorithm_exception(
        " HeatEquationParabolicPDESolverSystemParallelOMP::mult : An unknown operation mode was "
        "specified!");
  }
}

sgpp::base::DataVector* HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS() {
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

    applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
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
  } else {
    throw sgpp::base::algorithm_exception(
        "HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS : An unknown operation mode "
        "was specified!");
  }

  this->startTimestep();

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
  } else {
    throw sgpp::base::algorithm_exception(
        "HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS : An unknown operation mode "
        "was specified!");
  }

  rhs_complete.sub(result_complete);

  if (this->rhs != nullptr) {
    delete this->rhs;
  }

  this->rhs = new sgpp::base::DataVector(this->alpha_inner->getSize());
  this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

  return this->rhs;
}
}  // namespace pde
}  // namespace sgpp
