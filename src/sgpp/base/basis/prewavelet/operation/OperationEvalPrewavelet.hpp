/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de)
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)


#ifndef OPERATIONEVALPREWAVELET_HPP
#define OPERATIONEVALPREWAVELET_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
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

        virtual double eval(DataVector& alpha, std::vector<double>& point);
        virtual double test(DataVector& alpha, DataVector& data, DataVector& classes);
        virtual double integrate(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;

    };

  }
}

#endif /* OPERATIONEVELMODLINEAR_HPP */
