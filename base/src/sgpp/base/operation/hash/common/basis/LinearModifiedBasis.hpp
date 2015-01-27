// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODIFIED_LINEAR_BASE_HPP
#define MODIFIED_LINEAR_BASE_HPP

#include <cmath>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * modified linear base functions.
     *
     * @version $HEAD$
     */
    template<class LT, class IT>
    class LinearModifiedBasis: public Basis<LT, IT> {
      public:
        /**
         * Evaluate a basis function.
         * Has a dependence on the absolute position of grid point and support.
         */
        double eval(LT level, IT index, double p) {
          if (level == 1) {
            return 1.0;
          } else if (index == 1) {
            return 2.0 - (1 << level) * p;
          } else if (static_cast<int>(index) == static_cast<int>((1 << level) - 1)) {
            return (1 << level) * p - index + 1.0;
          }

          return 1.0 - fabs((1 << level) * p - index);
        }
    };

    // default type-def (unsigned int for level and index)
    typedef LinearModifiedBasis<unsigned int, unsigned int> SLinearModifiedBase;
  }
}

#endif /* MODIFIED_LINEAR_BASE_HPP */
