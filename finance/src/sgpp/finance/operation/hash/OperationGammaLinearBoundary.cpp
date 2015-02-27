// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationGammaLinearBoundary::OperationGammaLinearBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix& coef) : SGPP::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationGammaLinearBoundary::~OperationGammaLinearBoundary() {
    }

    void OperationGammaLinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLinearBoundary::upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * phi * dphi
      XPhidPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XPhidPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLinearBoundary::downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * phi * dphi
      XPhidPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XPhidPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLinearBoundary::upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLinearBoundary::downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLinearBoundary::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      SqXdPhidPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SqXdPhidPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLinearBoundary::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      SqXdPhidPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SqXdPhidPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}