// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonXLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationHestonXLinear::OperationHestonXLinear(sgpp::base::GridStorage* storage,
                                               sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationHestonXLinear::~OperationHestonXLinear() {}

void OperationHestonXLinear::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
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

void OperationHestonXLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                  size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::upOpDimOne(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::downOpDimOne(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  sgpp::finance::XdPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::finance::XdPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::upOpDimTwo(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::downOpDimTwo(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {}

void OperationHestonXLinear::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result, size_t dim) {}
}  // namespace finance
}  // namespace sgpp
