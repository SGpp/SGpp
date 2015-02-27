// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonELinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationHestonELinearBoundary::OperationHestonELinearBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef) : SGPP::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationHestonELinearBoundary::~OperationHestonELinearBoundary() {
    }

    void OperationHestonELinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);

    }

    void OperationHestonELinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);

    }

    void OperationHestonELinearBoundary::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonELinearBoundary::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}