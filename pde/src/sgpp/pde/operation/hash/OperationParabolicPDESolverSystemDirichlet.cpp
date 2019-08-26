// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationParabolicPDESolverSystemDirichlet.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationParabolicPDESolverSystemDirichlet::OperationParabolicPDESolverSystemDirichlet() {
  this->numSumGridpointsInner = 0;
  this->numSumGridpointsComplete = 0;
}

OperationParabolicPDESolverSystemDirichlet::~OperationParabolicPDESolverSystemDirichlet() {}

void OperationParabolicPDESolverSystemDirichlet::mult(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result) {
  result.setAll(0.0);

  if (this->tOperationMode == "ExEul") {
    applyMassMatrixInner(alpha, result);
  } else if (this->tOperationMode == "ImEul") {
    result.setAll(0.0);

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
    result.setAll(0.0);

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
  } else if (this->tOperationMode == "AdBas" || this->tOperationMode == "AdBasC") {
    result.setAll(0.0);

    applyMassMatrixInner(alpha, result);
  } else if (this->tOperationMode == "MPR") {
    applyMassMatrixInner(alpha, result);
  } else if (this->tOperationMode == "BDF2") {
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = (2.0 * tDiff + 1.0) / (tDiff + 1.0);
    result.setAll(0.0);

    sgpp::base::DataVector temp(alpha.getSize());

    applyMassMatrixInner(alpha, temp);

    temp.mult(alpha0);
    result.add(temp);

    applyLOperatorInner(alpha, temp);
    result.axpy((-1.0) * this->TimestepSize, temp);
  } else if (this->tOperationMode == "F23") {
    result.setAll(0.0);
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = 1.0 / (1.0 + tDiff);

    applyMassMatrixInner(alpha, result);
    result.mult(alpha0);

  } else {
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystem::mult : An unknown operation mode was specified!");
  }
}

sgpp::base::DataVector* OperationParabolicPDESolverSystemDirichlet::generateRHS() {
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
  } else if (this->tOperationMode == "AdBasC") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());

    applyMassMatrixComplete(*this->alpha_complete, temp);
    rhs_complete.add(temp);

    applyLOperatorComplete(*this->alpha_complete, temp);

    temp.mult((2.0) + this->TimestepSize / this->TimestepSize_old);

    sgpp::base::DataVector temp_old(this->alpha_complete->getSize());
    sgpp::base::DataVector alpha_old(this->alpha_complete->getSize());

    sgpp::base::DataVector temp_tmp(this->alpha_complete_old->getSize());
    temp_tmp.setAll(0.0);
    temp_tmp.add(*this->alpha_complete_old);

    double* OldData = alpha_old.getPointer();
    double* DataTmp = temp_tmp.getPointer();
    sgpp::base::GridStorage* gs = getGridStorage();
    sgpp::base::GridStorage* ogs = getOldGridStorage();
    sgpp::base::GridStorage::grid_map_iterator q;
    int length = 0;

    for (sgpp::base::GridStorage::grid_map_iterator p = gs->begin(); p != gs->end(); ++p) {
      q = ogs->find(p->first);

      if ((q->first)->equals(*p->first)) {
        size_t i = p->second;
        size_t j = q->second;

        if (j < temp_tmp.getSize()) {
          OldData[i] = DataTmp[j];
        } else {
          OldData[i] = 0;
        }
      } else {
        size_t i = p->second;
        OldData[i] = 0.0;
      }

      length++;
    }

    applyMassMatrixComplete(alpha_old, temp_old);

    applyLOperatorComplete(alpha_old, temp_old);

    temp_old.mult(this->TimestepSize / this->TimestepSize_old);

    temp.sub(temp_old);

    /*  DataVector temp3(this->alpha_complete->getSize());
      std::cout << "pre asd " << out_count++ << std::endl;

        std::cout << "AdBas size" << std::endl;
        std::cout << "AdBas size" << std::endl;

        double *OldData = temp_old.getPointer();
        double *Data = temp.getPointer();
        double *Data3 = temp3.getPointer();
        GridStorage *gs = getGridStorage();
        GridStorage *ogs = getOldGridStorage();
        GridStorage::grid_map_iterator q;
        int length = 0;
        for(GridStorage::grid_map_iterator p=gs->begin();p != gs->end();++p) {
          q = ogs->find(p->first);
          if((q->first)->equals(*p->first)) {
            int i = p->second;
            int j = q->second;
            Data3[length] = Data[i]-OldData[j];
            length++;

          }

        }
        temp3.resize(length);
        temp.resize(length);
        temp.setAll(0.0);
        temp.add(temp3);
    */

    rhs_complete.axpy((0.5) * this->TimestepSize, temp);
  } else if (this->tOperationMode == "MPR") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());

    sgpp::base::DataVector temp2(this->alpha_complete->getSize());

    applyMassMatrixComplete(*this->alpha_complete, temp);
    temp2.add(temp);

    applyLOperatorComplete(*this->alpha_complete, temp);
    temp2.axpy((0.5) * this->TimestepSize, temp);
    // rhs_complete.setAll(0.0);

    applyMassMatrixComplete(*this->alpha_complete, temp);
    rhs_complete.add(temp);

    applyLOperatorComplete(temp2, temp);
    rhs_complete.axpy((1.0) * this->TimestepSize, temp);

    //    applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
  } else if (this->tOperationMode == "BDF2") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());

    applyMassMatrixComplete(*this->alpha_complete, temp);

    double tDiff = this->TimestepSize / this->TimestepSize_old;

    double alpha1 = tDiff + 1.0;
    temp.mult(alpha1);
    rhs_complete.add(temp);

    sgpp::base::DataVector temp_old(this->alpha_complete->getSize());
    applyMassMatrixComplete(*this->alpha_complete_old, temp_old);

    double alpha2 = tDiff * tDiff / (1.0 + tDiff);
    temp_old.mult(alpha2);
    rhs_complete.sub(temp_old);
  } else if (this->tOperationMode == "F23") {
    rhs_complete.setAll(0.0);
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = (1.0 + tDiff);
    double alpha1 = alpha0 * (tDiff - 1.0);
    double alpha2 = -alpha0 * (tDiff * tDiff / (tDiff + 1.0));

    sgpp::base::DataVector temp(this->alpha_complete->getSize());
    sgpp::base::DataVector temp_old(this->alpha_complete->getSize());

    applyMassMatrixComplete(*this->alpha_complete, temp);
    temp.mult(alpha1);
    rhs_complete.sub(temp);

    applyLOperatorComplete(*this->alpha_complete, temp);
    temp.mult(alpha0 * this->TimestepSize);

    applyMassMatrixComplete(*this->alpha_complete_old, temp_old);
    temp_old.mult(alpha2);
    rhs_complete.sub(temp_old);

    rhs_complete.add(temp);
  } else {
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was "
        "specified!");
  }

  // Now we have the right hand side, lets apply the riskfree rate for the next timestep
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
  } else if (this->tOperationMode == "AdBas" || this->tOperationMode == "AdBasC") {
    applyMassMatrixComplete(alpha_bound, result_complete);
  } else if (this->tOperationMode == "MPR") {
    applyMassMatrixComplete(alpha_bound, result_complete);
  } else if (this->tOperationMode == "BDF2") {
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = (2.0 * tDiff + 1.0) / (tDiff + 1.0);
    sgpp::base::DataVector temp(alpha_bound.getSize());
    applyMassMatrixComplete(alpha_bound, temp);
    temp.mult(alpha0);
    result_complete.add(temp);

    applyLOperatorComplete(alpha_bound, temp);
    result_complete.axpy((-1.0) * this->TimestepSize, temp);

  } else if (this->tOperationMode == "F23") {
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = 1.0 / (1.0 + tDiff);

    applyMassMatrixComplete(alpha_bound, result_complete);
    result_complete.mult(alpha0);

  } else {
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was "
        "specified!");
  }

  rhs_complete.sub(result_complete);

  if (this->rhs != nullptr) {
    delete this->rhs;
  }

  this->rhs = new sgpp::base::DataVector(this->alpha_inner->getSize());
  this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

  return this->rhs;
}

sgpp::base::DataVector* OperationParabolicPDESolverSystemDirichlet::getGridCoefficientsForCG() {
  this->GridConverter->calcInnerCoefs(*this->alpha_complete, *this->alpha_inner);

  return this->alpha_inner;
}
}  // namespace pde
}  // namespace sgpp
