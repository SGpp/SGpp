/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/finance/basis/linearstretched/noboundary/operation/OperationDeltaLogLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiUpBBLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationDeltaLogLinearStretched::OperationDeltaLogLinearStretched(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef) : SGPP::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationDeltaLogLinearStretched::~OperationDeltaLogLinearStretched() {
    }

    void OperationDeltaLogLinearStretched::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinearStretched func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationDeltaLogLinearStretched::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinearStretched func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationDeltaLogLinearStretched::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearStretched func(this->storage);
      SGPP::base::sweep<DPhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationDeltaLogLinearStretched::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearStretched func(this->storage);
      SGPP::base::sweep<DPhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
