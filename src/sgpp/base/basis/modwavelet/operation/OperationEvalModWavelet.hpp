/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@mytum.de)

#ifndef OPERATIONEVALMODWAVELET_HPP
#define OPERATIONEVALMODWAVELET_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * This class implements OperationEval for a grids with mod wavelet basis ansatzfunctions with
     */
    class OperationEvalModWavelet : public OperationEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         */
        OperationEvalModWavelet(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationEvalModWavelet() {}

        virtual double eval(DataVector& alpha, std::vector<double>& point);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIINEVALMODWAVELET_HPP */
