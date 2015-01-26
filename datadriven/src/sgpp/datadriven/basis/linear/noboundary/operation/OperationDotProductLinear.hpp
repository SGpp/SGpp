/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONDOTPRODUCTLINEAR_HPP
#define OPERATIONDOTPRODUCTLINEAR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    class OperationDotProductLinear {
      public:
        /**
         * Constructor of OperationTestLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
    	OperationDotProductLinear(SGPP::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationDotProductLinear() {}

        virtual double eval(SGPP::base::DataVector& x1, SGPP::base::DataVector& x2);

      protected:
        /// Pointer to the grid's GridStorage object
        SGPP::base::GridStorage* storage;
    };

  }
}

#endif /* OPERATIONDOTPRODUCTLINEAR_HPP */
