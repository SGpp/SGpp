/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include <sgpp/pde/basis/linearstretched/noboundary/operation/OperationLaplaceLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/DowndPhidPhiBBIterativeLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationLaplaceLinearStretched::OperationLaplaceLinearStretched(SGPP::base::GridStorage* storage) : UpDownOneOpDim(storage) {
    }

    OperationLaplaceLinearStretched::~OperationLaplaceLinearStretched() {
    }

    void OperationLaplaceLinearStretched::specialOP(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim, size_t gradient_dim) {
      // In direction gradient_dim we only calculate the norm of the gradient
      // The up-part is empty, thus omitted
      if (dim > 0) {
        SGPP::base::DataVector temp(alpha.getSize());
        updown(alpha, temp, dim - 1, gradient_dim);
        downOpDim(temp, result, gradient_dim);
      } else {
        // Terminates dimension recursion
        downOpDim(alpha, result, gradient_dim);
      }
    }

    void OperationLaplaceLinearStretched::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      PhiPhiUpBBLinearStretched func(this->storage);
      SGPP::base::sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceLinearStretched::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      PhiPhiDownBBLinearStretched func(this->storage);
      SGPP::base::sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceLinearStretched::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      DowndPhidPhiBBIterativeLinearStretched myDown(this->storage);
      myDown(alpha, result, dim);
    }

    void OperationLaplaceLinearStretched::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
    }

  }
}
