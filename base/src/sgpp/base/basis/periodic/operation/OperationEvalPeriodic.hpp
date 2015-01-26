// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPERIODIC_HPP
#define OPERATIONEVALPERIODIC_HPP

#include <sgpp/base/operation/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements OperationEval for a grids with periodic linear basis ansatzfunctions with
     *
     * @version $HEAD$
     */
    class OperationEvalPeriodic : public OperationEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         */
    	OperationEvalPeriodic(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationEvalPeriodic() {}

        virtual double eval(DataVector& alpha, std::vector<double>& point);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;

    };

  }
}

#endif /* OPERATIONEVALPERIODIC_HPP */