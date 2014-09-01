/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_BASIS_BSPLINE_CLENSHAWCURTIS_OPERATION_OPERATIONNAIVEEVALBSPLINECLENSHAWCURTIS_HPP
#define SGPP_BASE_BASIS_BSPLINE_CLENSHAWCURTIS_OPERATION_OPERATIONNAIVEEVALBSPLINECLENSHAWCURTIS_HPP

#include "base/operation/OperationNaiveEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    /**
     * Operation for evaluating B-spline linear combinations on Clenshaw-Curtis grids.
     */
    class OperationNaiveEvalBsplineClenshawCurtis : public OperationNaiveEval {
      public:
        /**
         * Constructor.
         *
         * @param storage       storage of the sparse grid
         * @param degree        B-spline degree
         * @param cosine_table  cosine table for faster cosine evaluation (optional)
         */
        OperationNaiveEvalBsplineClenshawCurtis(GridStorage* storage, size_t degree,
                                                const CosineTable* cosine_table = NULL) :
          storage(storage),
          base(degree, cosine_table) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalBsplineClenshawCurtis() {
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
        /// 1D B-spline basis
        SBsplineClenshawCurtisBase base;
    };

  }
}

#endif
