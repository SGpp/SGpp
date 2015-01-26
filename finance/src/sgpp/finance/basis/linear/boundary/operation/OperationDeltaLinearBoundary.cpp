/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/finance/basis/linear/boundary/operation/OperationDeltaLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationDeltaLinearBoundary::OperationDeltaLinearBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef) : SGPP::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationDeltaLinearBoundary::~OperationDeltaLinearBoundary() {
    }

    void OperationDeltaLinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationDeltaLinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationDeltaLinearBoundary::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationDeltaLinearBoundary::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
