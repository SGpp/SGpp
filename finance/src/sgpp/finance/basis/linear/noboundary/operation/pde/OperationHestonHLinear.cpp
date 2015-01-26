/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#include <sgpp/finance/basis/linear/noboundary/operation/pde/OperationHestonHLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationHestonHLinear::OperationHestonHLinear(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix& coef) : SGPP::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationHestonHLinear::~OperationHestonHLinear() {
    }

    void OperationHestonHLinear::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
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

    void OperationHestonHLinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonHLinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonHLinear::upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonHLinear::downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonHLinear::upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * phi * phi
      XPhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<XPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonHLinear::downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * phi * phi
      XPhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<XPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonHLinear::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
    }

    void OperationHestonHLinear::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
    }

  }
}
