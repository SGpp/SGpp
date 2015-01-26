// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/boundary/operation/OperationLFLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationLFLinearBoundary::OperationLFLinearBoundary(SGPP::base::GridStorage* storage) : SGPP::pde::StdUpDown(storage) {
    }

    OperationLFLinearBoundary::~OperationLFLinearBoundary() {
    }

    void OperationLFLinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // X * dphi * phi
      XdPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLFLinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // X * dphi * phi
      XdPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}