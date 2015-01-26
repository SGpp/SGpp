/* ****************************************************************************
* Copyright (C) 2008-2010 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de)
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONCONVERT_HPP
#define OPERATIONCONVERT_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#ifdef _WIN32
#pragma warning(disable: 4267)
#endif

namespace sg {
  namespace base {

    /**
     * Operation that cpnverts a given basis into the normal, linear hat basis and vice versa
     *
     */
    class OperationConvert {
      public:
        /**
         * Constructor
         */
        OperationConvert() {}

        /**
         * Destructor
         */
        virtual ~OperationConvert() {}

        /**
         * Convert given basis into linear hat basis.
         */
        virtual void doConvertToLinear(DataVector& alpha) = 0;


        /**
         * Convert from a linear coefficient vector into given basis.
         */
        virtual void doConvertFromLinear(DataVector& alpha) = 0;
    };

  }
}

#endif /* OPERATIONCONVERT_HPP */
