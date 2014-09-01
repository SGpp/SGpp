/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_WAVELET_BOUNDARY_OPERATION_OPERATIONNAIVEEVALWAVELETBOUNDARY_HPP
#define SGPP_BASE_BASIS_WAVELET_BOUNDARY_OPERATION_OPERATIONNAIVEEVALWAVELETBOUNDARY_HPP

#include "base/operation/OperationNaiveEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    /**
     * Operation for evaluating wavelet linear combinations on Boundary grids.
     */
    class OperationNaiveEvalWaveletBoundary : public OperationNaiveEval {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         */
        OperationNaiveEvalWaveletBoundary(GridStorage* storage) : storage(storage) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalWaveletBoundary() {
        }

        /**
         * @param alpha     coefficient vector
         * @param point     evaluation point
         * @return          value of linear combination
         */
        virtual double eval(DataVector& alpha, std::vector<double>& point);

      protected:
        /// storage of the sparse grid
        GridStorage* storage;
        /// 1D wavelet basis
        SWaveletBoundaryBase base;
    };

  }
}

#endif
