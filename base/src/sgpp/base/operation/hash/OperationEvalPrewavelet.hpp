// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPREWAVELET_HPP
#define OPERATIONEVALPREWAVELET_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements OperationEval for a grids with prewavelet basis ansatzfunctions without boundaries
     */
    class OperationEvalPrewavelet : public OperationEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         */
        OperationEvalPrewavelet(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationEvalPrewavelet() {}

        virtual float_t eval(DataVector& alpha, std::vector<float_t>& point);
        virtual float_t test(DataVector& alpha, DataVector& data, DataVector& classes);
        virtual float_t integrate(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;

    };

  }
}

#endif /* OPERATIONEVELMODLINEAR_HPP */