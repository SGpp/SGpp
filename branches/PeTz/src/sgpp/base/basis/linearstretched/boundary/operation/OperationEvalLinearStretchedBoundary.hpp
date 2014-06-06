/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef OPERATIONEVALLINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONEVALLINEARSTRETCHEDBOUNDARY_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * This class implements OperationEval for a grids with linear basis ansatzfunctions with
     * boundaries
     *
     * @version $HEAD$
     */
    class OperationEvalLinearStretchedBoundary : public OperationEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         */
        OperationEvalLinearStretchedBoundary(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationEvalLinearStretchedBoundary() {}

        virtual double eval(DataVector& alpha, std::vector<double>& point);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONEVALLINEARSTRETCHEDBOUNDARY_HPP */
