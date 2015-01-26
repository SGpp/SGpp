/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef OPERATIONEVALLINEARSTRETCHED_HPP
#define OPERATIONEVALLINEARSTRETCHED_HPP

#include <sgpp/base/operation/OperationEval.hpp>
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
        virtual ~OperationEvalLinearStretched() {}

        virtual double eval(DataVector& alpha, std::vector<double>& point);

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONEVALLINEARSTRETCHED_HPP */
