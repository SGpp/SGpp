// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonYLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiUpBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiDownBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationHestonYLinearBoundary::OperationHestonYLinearBoundary(sgpp::base::GridStorage* storage,
                                                               sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationHestonYLinearBoundary::~OperationHestonYLinearBoundary() {}

void OperationHestonYLinearBoundary::mult(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result) {
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

void OperationHestonYLinearBoundary::up(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonYLinearBoundary::down(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonYLinearBoundary::upOpDimOne(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  sgpp::finance::SqXdPhidPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::finance::SqXdPhidPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonYLinearBoundary::downOpDimOne(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  sgpp::finance::SqXdPhidPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::finance::SqXdPhidPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonYLinearBoundary::upOpDimTwo(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonYLinearBoundary::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // x * phi * phi
  XPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonYLinearBoundary::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                           sgpp::base::DataVector& result,
                                                           size_t dim) {}

void OperationHestonYLinearBoundary::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result,
                                                             size_t dim) {}
}  // namespace finance
}  // namespace sgpp
