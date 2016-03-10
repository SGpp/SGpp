// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonWLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationHestonWLinearBoundary::OperationHestonWLinearBoundary(sgpp::base::GridStorage* storage,
                                                               sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationHestonWLinearBoundary::~OperationHestonWLinearBoundary() {}

void OperationHestonWLinearBoundary::mult(sgpp::base::DataVector& alpha,
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

void OperationHestonWLinearBoundary::up(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonWLinearBoundary::down(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonWLinearBoundary::upOpDimOne(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // X * phi * dphi
  XPhidPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhidPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonWLinearBoundary::downOpDimOne(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // X * phi * dphi
  XPhidPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhidPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonWLinearBoundary::upOpDimTwo(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // X * dphi * phi
  XdPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonWLinearBoundary::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // X * dphi * phi
  XdPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonWLinearBoundary::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                           sgpp::base::DataVector& result,
                                                           size_t dim) {}

void OperationHestonWLinearBoundary::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result,
                                                             size_t dim) {}
}  // namespace finance
}  // namespace sgpp
