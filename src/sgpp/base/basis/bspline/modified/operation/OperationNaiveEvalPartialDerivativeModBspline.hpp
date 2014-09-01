/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_BSPLINE_MODIFIED_OPERATION_OPERATIONEVALPARTIALDERIVATIVEMODBSPLINE_HPP
#define SGPP_BASE_BASIS_BSPLINE_MODIFIED_OPERATION_OPERATIONEVALPARTIALDERIVATIVEMODBSPLINE_HPP

#include "base/operation/OperationNaiveEvalPartialDerivative.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/modified/ModBsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    /**
     * Operation for evaluating partial derivatives of modified B-spline
     * linear combinations on Noboundary grids.
     */
    class OperationNaiveEvalPartialDerivativeModBspline :
      public OperationNaiveEvalPartialDerivative {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         * @param degree    B-spline degree
         */
        OperationNaiveEvalPartialDerivativeModBspline(GridStorage* storage, size_t degree) :
          storage(storage), base(degree) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalPartialDerivativeModBspline() {
        }

        /**
         * @param alpha     coefficient vector
         * @param point     evaluation point
         * @param deriv_dim dimension in which the partial derivative should be taken
         * @return          value of the partial derivative of the linear combination
         */
        virtual double evalPartialDerivative(DataVector& alpha,
                                             const std::vector<double>& point,
                                             size_t deriv_dim);

      protected:
        /// storage of the sparse grid
        GridStorage* storage;
        /// 1D B-spline basis
        SModBsplineBase base;
    };

  }
}

#endif
