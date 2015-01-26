/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include <sgpp/finance/basis/linearstretched/noboundary/operation/OperationDeltaLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiUpBBLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationDeltaLinearStretched::OperationDeltaLinearStretched(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef) : SGPP::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationDeltaLinearStretched::~OperationDeltaLinearStretched() {
    }

    void OperationDeltaLinearStretched::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinearStretched func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationDeltaLinearStretched::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinearStretched func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationDeltaLinearStretched::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiUpBBLinearStretched func(this->storage);
      SGPP::base::sweep<XdPhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationDeltaLinearStretched::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiDownBBLinearStretched func(this->storage);
      SGPP::base::sweep<XdPhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
