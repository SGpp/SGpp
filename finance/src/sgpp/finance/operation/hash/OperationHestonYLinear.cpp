// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonYLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationHestonYLinear::OperationHestonYLinear(sgpp::base::GridStorage* storage,
                                               sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationHestonYLinear::~OperationHestonYLinear() {}

void OperationHestonYLinear::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
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

void OperationHestonYLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                  size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::upOpDimOne(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  sgpp::finance::SqXdPhidPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::finance::SqXdPhidPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::downOpDimOne(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  sgpp::finance::SqXdPhidPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::finance::SqXdPhidPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::upOpDimTwo(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::downOpDimTwo(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {}

void OperationHestonYLinear::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result, size_t dim) {}
}  // namespace finance
}  // namespace sgpp
