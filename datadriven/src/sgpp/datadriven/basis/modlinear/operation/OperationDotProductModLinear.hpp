/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONDOTPRODUCTMODLINEAR_HPP
#define OPERATIONDOTPRODUCTMODLINEAR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sg {
  namespace datadriven {

    class OperationDotProductModLinear {
      public:
        /**
         * Constructor of OperationTestLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
    	OperationDotProductModLinear(sg::base::GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationDotProductModLinear() {}

        virtual double eval(sg::base::DataVector& x1, sg::base::DataVector& x2);

      protected:
        /// Pointer to the grid's GridStorage object
        sg::base::GridStorage* storage;
    };

  }
}

#endif /* OPERATIONDOTPRODUCTMODLINEAR_HPP */
