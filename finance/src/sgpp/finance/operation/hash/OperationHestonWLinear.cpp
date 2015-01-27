// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonWLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationHestonWLinear::OperationHestonWLinear(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix& coef) : SGPP::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationHestonWLinear::~OperationHestonWLinear() {
    }

    void OperationHestonWLinear::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
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
                    {
                      result.axpy(this->coefs->get(i, j), beta);
                    }
                  }
                } else {
                  this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j);

                  #pragma omp critical
                  {
                    result.add(beta);
                  }
                }
              }
            }
          }

          #pragma omp taskwait
        }
      }
    }

    void OperationHestonWLinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonWLinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonWLinear::upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // X * phi * dphi
      XPhidPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<XPhidPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonWLinear::downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // X * phi * dphi
      XPhidPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<XPhidPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonWLinear::upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // X * dphi * phi
      XdPhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonWLinear::downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // X * dphi * phi
      XdPhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonWLinear::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
    }

    void OperationHestonWLinear::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
    }

  }
}