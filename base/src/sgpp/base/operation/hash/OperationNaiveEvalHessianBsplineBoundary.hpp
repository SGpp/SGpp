// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALHESSIANBSPLINEBOUNDARY_HPP
#define OPERATIONEVALHESSIANBSPLINEBOUNDARY_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessian.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace SGPP {
  namespace base {

    /**
     * Operation for evaluating B-spline linear combinations on Boundary grids, their gradients
     * and their Hessians.
     */
    class OperationNaiveEvalHessianBsplineBoundary : public OperationNaiveEvalHessian {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         * @param degree    B-spline degree
         */
        OperationNaiveEvalHessianBsplineBoundary(GridStorage* storage, size_t degree) :
          storage(storage), base(degree) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalHessianBsplineBoundary() {
        }

        /**
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient vector of linear combination
         * @param[out]  hessian     Hessian matrix of linear combination
         * @return                  value of linear combination
         */
        virtual float_t evalHessian(DataVector& alpha, const std::vector<float_t>& point,
                                    DataVector& gradient, DataMatrix& hessian);

      protected:
        /// storage of the sparse grid
        GridStorage* storage;
        /// 1D B-spline basis
        SBsplineBoundaryBase base;
    };

  }
}

#endif /* OPERATIONEVALHESSIANBSPLINEBOUNDARY_HPP */
