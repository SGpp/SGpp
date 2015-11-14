// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALGRADIENTWAVELETBOUNDARY_HPP
#define OPERATIONNAIVEEVALGRADIENTWAVELETBOUNDARY_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace base {

    /**
     * Operation for evaluating wavelet linear combinations on Boundary grids and their gradients.
     */
    class OperationNaiveEvalGradientWaveletBoundary : public OperationNaiveEvalGradient {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         */
        OperationNaiveEvalGradientWaveletBoundary(GridStorage* storage) : storage(storage) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalGradientWaveletBoundary() {
        }

        /**
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient of linear combination
         * @return                  value of linear combination
         */
        virtual float_t evalGradient(const DataVector& alpha,
                                     const DataVector& point,
                                     DataVector& gradient);

      protected:
        /// storage of the sparse grid
        GridStorage* storage;
        /// 1D wavelet basis
        SWaveletBoundaryBase base;
    };

  }
}

#endif /* OPERATIONNAIVEEVALGRADIENTWAVELETBOUNDARY_HPP */
