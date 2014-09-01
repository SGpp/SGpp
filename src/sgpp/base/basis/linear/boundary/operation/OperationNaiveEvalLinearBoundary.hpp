/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_LINEAR_BOUNDARY_OPERATION_OPERATIONNAIVEEVALLINEARBOUNDARY_HPP
#define SGPP_BASE_BASIS_LINEAR_BOUNDARY_OPERATION_OPERATIONNAIVEEVALLINEARBOUNDARY_HPP

#include "base/operation/OperationNaiveEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    class OperationNaiveEvalLinearBoundary : public OperationNaiveEval {
      public:
        /**
         * Constructor.
         *
         * @param storage   storage of the sparse grid
         */
        OperationNaiveEvalLinearBoundary(GridStorage* storage) : storage(storage) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalLinearBoundary() {
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
        /// 1D linear basis
        SLinearBoundaryBase base;
    };

  }
}

#endif
