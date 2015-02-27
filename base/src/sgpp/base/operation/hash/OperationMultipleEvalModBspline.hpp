// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALMODBSPLINE_HPP
#define OPERATIONMULTIPLEEVALMODBSPLINE_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements OperationMultipleEval for a grid with modified Bspline basis functions
     *
     * @version $HEAD$
     */
    class OperationMultipleEvalModBspline: public OperationMultipleEval {
      public:
        /**
         * Constructor
         *
         * @param grid grid
         * @param degree the Bspline's degree
         * @param dataset the dataset that should be evaluated
         */
        OperationMultipleEvalModBspline(Grid& grid, size_t degree,
                                        DataMatrix& dataset) :
          OperationMultipleEval(grid, dataset), base(degree) {
          this->storage = grid.getStorage();
        }

        /**
         * Destructor
         */
        virtual ~OperationMultipleEvalModBspline() {
        }

        virtual void mult(DataVector& alpha, DataVector& result);
        virtual void multTranspose(DataVector& source, DataVector& result);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
        /// Mod Bspline Basis object
        SBsplineModifiedBase base;
    };

  }
}

#endif /* OPERATIONMULTIPLEEVALMODBSPLINE_HPP */
