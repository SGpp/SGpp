// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPARTIALDERIVATIVEMODBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONEVALPARTIALDERIVATIVEMODBSPLINECLENSHAWCURTIS_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivative.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace base {

    /**
     * Operation for evaluating partial derivatives of
     * modified Clenshaw-Curtis B-spline
     * linear combinations on Noboundary grids.
     */
    class OperationNaiveEvalPartialDerivativeModBsplineClenshawCurtis :
      public OperationNaiveEvalPartialDerivative {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         * @param degree    B-spline degree
         */
        OperationNaiveEvalPartialDerivativeModBsplineClenshawCurtis(
          GridStorage* storage, size_t degree) :
          storage(storage), base(degree) {
        }

        /**
         * Destructor.
         */
        virtual ~OperationNaiveEvalPartialDerivativeModBsplineClenshawCurtis() override {
        }

        /**
         * @param alpha     coefficient vector
         * @param point     evaluation point
         * @param derivDim  dimension in which the partial derivative should be taken
         * @return          value of the partial derivative of the linear combination
         */
        virtual float_t evalPartialDerivative(const DataVector& alpha,
                                              const DataVector& point,
                                              size_t derivDim) override;

      protected:
        /// storage of the sparse grid
        GridStorage* storage;
        /// 1D B-spline basis
        SBsplineModifiedClenshawCurtisBase base;
    };

  }
}

#endif /* OPERATIONEVALPARTIALDERIVATIVEMODBSPLINECLENSHAWCURTIS_HPP */
