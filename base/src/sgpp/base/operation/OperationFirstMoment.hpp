/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)
// @author Benjamin

#ifndef OPERATIONFIRSTMOMENT_HPP
#define OPERATIONFIRSTMOMENT_HPP

#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    /**
     * This class provides the first moment of a sparse grid function
     */
    class OperationFirstMoment {
      public:
        /**
         * Constructor
         */
        OperationFirstMoment() {}

        /**
         * Destructor
         */
        virtual ~OperationFirstMoment() {}

        /**
         * Integrate the sparse grid function
         *
         * @param alpha the function's values in the nodal basis
         */
        virtual double doQuadrature(DataVector& alpha) = 0;

    };

  }
}

#endif /* OPERATIONFIRSTMOMENT_HPP */
