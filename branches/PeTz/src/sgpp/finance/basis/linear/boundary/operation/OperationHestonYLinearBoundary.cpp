/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/boundary/operation/OperationHestonYLinearBoundary.hpp"

#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/XPhiPhiUpBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/XPhiPhiDownBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp"


#include "base/algorithm/sweep.hpp"

#include <iostream>

namespace sg {
  namespace finance {

    OperationHestonYLinearBoundary::OperationHestonYLinearBoundary(sg::base::GridStorage* storage, sg::base::DataMatrix& coef) : sg::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationHestonYLinearBoundary::~OperationHestonYLinearBoundary() {
    }

    void OperationHestonYLinearBoundary::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
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
                sg::base::DataVector beta(result.getSize());

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

    void OperationHestonYLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonYLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonYLinearBoundary::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      sg::finance::SqXdPhidPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::finance::SqXdPhidPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonYLinearBoundary::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      sg::finance::SqXdPhidPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::finance::SqXdPhidPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonYLinearBoundary::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * phi * phi
      XPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<XPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonYLinearBoundary::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * phi * phi
      XPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<XPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonYLinearBoundary::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
    }

    void OperationHestonYLinearBoundary::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
    }

  }
}
