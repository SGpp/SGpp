/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_WAVELET_NOBOUNDARY_OPERATION_OPERATIONNAIVEEVALGRADIENTWAVELET_HPP
#define SGPP_BASE_BASIS_WAVELET_NOBOUNDARY_OPERATION_OPERATIONNAIVEEVALGRADIENTWAVELET_HPP

#include "base/operation/OperationNaiveEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    /**
     * Operation for evaluating wavelet linear combinations on Noboundary grids and their gradients.
     */
    class OperationNaiveEvalGradientWavelet : public OperationNaiveEvalGradient {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         */
        OperationNaiveEvalGradientWavelet(GridStorage* storage) : storage(storage) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalGradientWavelet() {
        }

        /**
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient of linear combination
         * @return                  value of linear combination
         */
        virtual double evalGradient(DataVector& alpha, const std::vector<double>& point,
                                    DataVector& gradient);

      protected:
        /// storage of the sparse grid
        GridStorage* storage;
        /// 1D wavelet basis
        SWaveletBase base;
    };

  }
}

#endif
