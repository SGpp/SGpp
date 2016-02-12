// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonBLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/pde/operation/hash/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/UpdPhidPhiBBIterativeLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace SGPP {
namespace finance {

OperationHestonBLinearBoundary::OperationHestonBLinearBoundary(SGPP::base::GridStorage* storage,
                                                               SGPP::base::DataMatrix& coef)
    : SGPP::pde::UpDownTwoOpDims(storage, coef) {}

OperationHestonBLinearBoundary::~OperationHestonBLinearBoundary() {}

void OperationHestonBLinearBoundary::mult(SGPP::base::DataVector& alpha,
                                          SGPP::base::DataVector& result) {
  result.setAll(0.0);

#pragma omp parallel
  {
#pragma omp single nowait
    {
      for (size_t i = 0; i < this->numAlgoDims_; i++) {
        for (size_t j = 0; j < this->numAlgoDims_; j++) {
// no symmetry in the operator
#pragma omp task firstprivate(i, j) shared(alpha, result)
          {
            SGPP::base::DataVector beta(result.getSize());

            if (this->coefs != NULL) {
              if (this->coefs->get(i, j) != 0.0) {
                this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j);

#pragma omp critical
                { result.axpy(this->coefs->get(i, j), beta); }
              }
            } else {
              this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j);

#pragma omp critical
              { result.add(beta); }
            }
          }
        }
      }

#pragma omp taskwait
    }
  }
}

void OperationHestonBLinearBoundary::up(SGPP::base::DataVector& alpha,
                                        SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonBLinearBoundary::down(SGPP::base::DataVector& alpha,
                                          SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonBLinearBoundary::upOpDimOne(SGPP::base::DataVector& alpha,
                                                SGPP::base::DataVector& result, size_t dim) {
  SGPP::pde::UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
  myUp(alpha, result, dim);
}

void OperationHestonBLinearBoundary::downOpDimOne(SGPP::base::DataVector& alpha,
                                                  SGPP::base::DataVector& result, size_t dim) {
  // Dphi * dphi
  SGPP::pde::DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
  myDown(alpha, result, dim);
}

void OperationHestonBLinearBoundary::upOpDimTwo(SGPP::base::DataVector& alpha,
                                                SGPP::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiUpBBLinearBoundary func(this->storage);
  SGPP::base::sweep<XPhiPhiUpBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonBLinearBoundary::downOpDimTwo(SGPP::base::DataVector& alpha,
                                                  SGPP::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiDownBBLinearBoundary func(this->storage);
  SGPP::base::sweep<XPhiPhiDownBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonBLinearBoundary::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                           SGPP::base::DataVector& result,
                                                           size_t dim) {}

void OperationHestonBLinearBoundary::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                             SGPP::base::DataVector& result,
                                                             size_t dim) {}
}  // namespace finance
}  // namespace SGPP
