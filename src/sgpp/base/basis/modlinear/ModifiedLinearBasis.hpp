/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef MODIFIED_LINEAR_BASE_HPP
#define MODIFIED_LINEAR_BASE_HPP

#include "base/basis/Basis.hpp"

#include <cmath>

namespace sg {
  namespace base {

    /**
     * modified linear base functions.
     *
     * @version $HEAD$
     */
    template<class LT, class IT>
    class ModifiedLinearBasis : public Basis<LT, IT> {
      public:
        /**
         * Evaluate a basis function.
         * Has a dependence on the absolute position of grid point and support.
         */
        double eval(LT level, IT index, double p) {
          double hinv = static_cast<double>(1 << level);
          double h = 1.0 / hinv;
          
          if (level == 1) {
            return 1.0;
          } else if (index == 1) {
            return ((p <= 2.0 * h) ? (2.0 - hinv * p) : 0.0);
          } else if (static_cast<int>(index) == static_cast<int>((1 << level) - 1)) {
            return ((p >= 1.0 - 2.0 * h) ? (hinv * p - index + 1.0) : 0.0);
          }
          
          double result = 1.0 - fabs(hinv * p - index);
          
          if (result >= 0.0) {
            return result;
          } else {
            return 0.0;
          }
        }
    };

    // default type-def (unsigned int for level and index)
    typedef ModifiedLinearBasis<unsigned int, unsigned int> SModLinearBase;
  }
}

#endif /* MODIFIED_LINEAR_BASE_HPP */

