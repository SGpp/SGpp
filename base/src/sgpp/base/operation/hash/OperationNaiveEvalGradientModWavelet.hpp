// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALGRADIENTMODWAVELET_HPP
#define OPERATIONNAIVEEVALGRADIENTMODWAVELET_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace base {

    /**
     * Operation for evaluating modified wavelet linear combinations on Noboundary grids and
     * their gradients.
     */
    class OperationNaiveEvalGradientModWavelet : public OperationNaiveEvalGradient {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         */
        OperationNaiveEvalGradientModWavelet(GridStorage* storage) : storage(storage) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalGradientModWavelet() {
        }

        /**
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient of linear combination
         * @return                  value of linear combination
         */
        virtual float_t evalGradient(DataVector& alpha, const std::vector<float_t>& point,
                                    DataVector& gradient);

      protected:
        /// storage of the sparse grid
        GridStorage* storage;
        /// 1D wavelet basis
        SWaveletModifiedBase base;
    };

  }
}

#endif /* OPERATIONNAIVEEVALGRADIENTMODWAVELET_HPP */
