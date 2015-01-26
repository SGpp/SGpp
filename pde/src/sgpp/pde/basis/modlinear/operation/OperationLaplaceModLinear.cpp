/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/basis/modlinear/operation/OperationLaplaceModLinear.hpp>

#include <sgpp/pde/basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp>
#include <sgpp/pde/basis/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp>
#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp>
#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiUpModLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationLaplaceModLinear::OperationLaplaceModLinear(SGPP::base::GridStorage* storage) : UpDownOneOpDim(storage) {
    }

    OperationLaplaceModLinear::~OperationLaplaceModLinear() {
    }

    void OperationLaplaceModLinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      PhiPhiUpModLinear func(this->storage);
      SGPP::base::sweep<PhiPhiUpModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceModLinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      PhiPhiDownModLinear func(this->storage);
      SGPP::base::sweep<PhiPhiDownModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceModLinear::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      dPhidPhiDownModLinear func(this->storage);
      SGPP::base::sweep<dPhidPhiDownModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceModLinear::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      dPhidPhiUpModLinear func(this->storage);
      SGPP::base::sweep<dPhidPhiUpModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

  }
}

