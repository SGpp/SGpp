// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonHLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationHestonHLinear::OperationHestonHLinear(sgpp::base::GridStorage* storage,
                                               sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationHestonHLinear::~OperationHestonHLinear() {}

void OperationHestonHLinear::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
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
            sgpp::base::DataVector beta(result.getSize());

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

void OperationHestonHLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                  size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::upOpDimOne(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::downOpDimOne(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::upOpDimTwo(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::downOpDimTwo(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {}

void OperationHestonHLinear::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result, size_t dim) {}
}  // namespace finance
}  // namespace sgpp
