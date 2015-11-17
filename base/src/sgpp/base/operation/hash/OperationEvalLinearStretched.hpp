// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALLINEARSTRETCHED_HPP
#define OPERATIONEVALLINEARSTRETCHED_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements OperationEval for a grids with linear basis ansatzfunctions without boundaries
     */
    class OperationEvalLinearStretched : public OperationEval {
      public:
        /**
         * Constructor of OperationEvalLinearStretched
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationEvalLinearStretched(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationEvalLinearStretched() override {}

        virtual float_t eval(const DataVector& alpha,
                             const DataVector& point) override;

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONEVALLINEARSTRETCHED_HPP */
