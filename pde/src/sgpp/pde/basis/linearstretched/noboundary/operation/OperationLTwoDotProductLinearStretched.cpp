/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include <sgpp/pde/basis/linearstretched/noboundary/operation/OperationLTwoDotProductLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationLTwoDotProductLinearStretched::OperationLTwoDotProductLinearStretched(SGPP::base::GridStorage* storage) : StdUpDown(storage) {
    }

    OperationLTwoDotProductLinearStretched::~OperationLTwoDotProductLinearStretched() {
    }

    void OperationLTwoDotProductLinearStretched::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiUpBBLinearStretched func(this->storage);
      SGPP::base::sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationLTwoDotProductLinearStretched::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiDownBBLinearStretched func(this->storage);
      SGPP::base::sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
